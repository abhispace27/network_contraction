#!/usr/bin/env python3
"""
Force/stress analysis pipeline for circular elastic-network simulations.

This refactor keeps the original scientific workflow, but organizes it as a
maintainable Python script suitable for a GitHub portfolio:

1. Load node forces, positions, strains, rest lengths, connectivity, boundary nodes, and dipole nodes.
2. Reconstruct stretching forces from bond strain.
3. Estimate bending-force contributions.
4. Project boundary forces along radial directions.
5. Compute far-field and local force-dipole moments.
6. Compute radial stress on outer/inner boundaries and concentric rings.
7. Extract geometric features of dipole-center organization.
8. Save tabular outputs and summary plots.

Run example:
    python force_stress_analysis_refactored.py \
        --base-dir /home/abhinav/david \
        --lattice-size 64 \
        --num-centers 5 \
        --num-force-steps 10 \
        --pbond 0.55 \
        --kappa 1e-6 \
        --seed 112 \
        --network-id 111111,101111 \
        --run-type cluster_radial
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import ConvexHull

RunType = Literal[
    "default",
    "inner_fixed",
    "radial",
    "inner_radial",
    "inner_radial_arp",
    "radial_only",
    "cluster",
    "cluster_bash",
    "cluster_hex",
    "cluster_random_hex",
    "cluster_radial",
    "force",
]


# =============================================================================
# Configuration and filename handling
# =============================================================================


@dataclass(frozen=True)
class AnalysisConfig:
    base_dir: Path
    lattice_size: int = 64
    num_centers: int = 5
    num_force_steps: int = 10
    pbond: float = 0.55
    tolerance: float = 1.0e-7
    kappa: float = 1.0e-6
    rest_length: float = 0.9
    mu_tension: float = 1.0
    mu_compression: float = 1.0
    seed: int = 112
    network_id: str = "111111,101111"
    ranseed: str = "667720,601210"
    run_type: RunType = "cluster_radial"
    use_seed: bool = True
    only_last: bool = False
    only_last_num: int = 2
    force_value: float = 0.01
    dont_throw_disconnected_dipoles: bool = False
    center_node_64: int = 2080
    center_node_128: int = 8383
    inner_boundary_radius: float = 13.0
    inner_boundary_width: float = 0.5
    outer_radius: int = 25
    ring_width: float = 0.5

    @property
    def num_dipoles(self) -> int:
        return 6 * self.num_centers

    @property
    def unconnected_value(self) -> int:
        return 9999 if self.lattice_size == 64 else -1

    @property
    def center_node(self) -> int:
        return self.center_node_64 if self.lattice_size == 64 else self.center_node_128

    @property
    def pbond_str(self) -> str:
        return f"{self.pbond:.2f}"

    @property
    def tol_str(self) -> str:
        return f"{self.tolerance:.2e}"

    @property
    def kappa_str(self) -> str:
        return "0.00e+00" if self.kappa == 0 else f"{self.kappa:.2e}"

    @property
    def rest_length_str(self) -> str:
        return f"{self.force_value:.9f}" if self.run_type == "force" else f"{self.rest_length:.4f}"

    @property
    def mu_str(self) -> str:
        return f"{self.mu_tension:.4f}"

    @property
    def mu_c_str(self) -> str:
        return f"{self.mu_compression:.4f}"

    @property
    def final_step(self) -> int:
        return self.num_force_steps

    @property
    def position_steps(self) -> range:
        if self.only_last:
            return range(self.only_last_num)
        return range(self.num_force_steps + 1)

    @property
    def force_steps(self) -> range:
        if self.only_last:
            return range(1, self.only_last_num)
        return range(1, self.num_force_steps + 1)

    @property
    def rest_length_steps(self) -> np.ndarray:
        step = (1.0 - self.rest_length) / self.num_force_steps
        return np.array([1.0 - i * step for i in range(1, self.num_force_steps + 1)])


class PathBuilder:
    def __init__(self, cfg: AnalysisConfig):
        self.cfg = cfg
        self.folder = self._build_folder()
        self.centers_folder = self._build_centers_folder()

    def _kappa_folder(self) -> str:
        if self.cfg.kappa == 0:
            return "kappa2_0"
        if np.isclose(self.cfg.kappa, 1e-6):
            return "kappa2_e-6"
        if np.isclose(self.cfg.kappa, 1e-5):
            return "kappa2_e-5"
        if np.isclose(self.cfg.kappa, 1e-4):
            return "kappa2_e-4"
        return f"kappa2_{self.cfg.kappa:.0e}"

    def _root_for_run_type(self) -> str:
        mapping: dict[str, str] = {
            "default": "lattice_nopbd_bndry_all_clamp_restlength_circle",
            "inner_fixed": "lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed",
            "radial": "lattice_nopbd_bndry_all_clamp_restlength_circle_radial",
            "inner_radial": "lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial",
            "inner_radial_arp": "lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp",
            "radial_only": "lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_only_radial_arp",
            "cluster": "lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_cluster",
            "cluster_bash": "cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash",
            "cluster_hex": "cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex",
            "cluster_random_hex": "cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds",
            "cluster_radial": "cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial",
            "force": "cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_force",
        }
        return mapping[self.cfg.run_type]

    def _build_folder(self) -> Path:
        if self.cfg.pbond == 1:
            return Path(self._root_for_run_type())
        return (
            Path(self._root_for_run_type())
            / self._kappa_folder()
            / self.cfg.ranseed
            / str(self.cfg.seed)
            / self.cfg.network_id
        )

    def _build_centers_folder(self) -> Path:
        if self.cfg.pbond == 1:
            return self.folder
        return self.folder.parent

    def input_dir(self, subdir: str) -> Path:
        return self.cfg.base_dir / self.folder / "txt" / subdir

    def output_dir(self, subdir: str) -> Path:
        path = self.cfg.base_dir / self.folder / "txt" / subdir
        path.mkdir(parents=True, exist_ok=True)
        return path

    def means_dir(self) -> Path:
        path = self.cfg.base_dir / self.centers_folder / "means"
        path.mkdir(parents=True, exist_ok=True)
        return path

    def plot_dir(self, subdir: str = "strain_hist") -> Path:
        path = self.cfg.base_dir / self.folder / "png" / subdir
        path.mkdir(parents=True, exist_ok=True)
        return path

    def base_tag(self, *, include_double_underscore: bool = True) -> str:
        cfg = self.cfg
        pieces: list[str] = []

        if cfg.use_seed and (cfg.pbond == 1 or cfg.run_type == "force"):
            pieces.append(f"srand_{cfg.seed}")

        if cfg.lattice_size != 64:
            pieces.append(str(cfg.lattice_size))

        pieces.extend(
            [
                str(cfg.num_centers),
                str(cfg.num_dipoles),
                cfg.pbond_str,
                cfg.tol_str,
                cfg.kappa_str,
            ]
        )

        left = "_".join(pieces)
        right = "_".join([cfg.rest_length_str, cfg.mu_str, cfg.mu_c_str])

        if include_double_underscore:
            return f"{left}__{right}"
        return f"{left}_{right}"

    def stepped_tag(self, step: int, *, suffix_force: bool = True) -> str:
        tag = f"{self.base_tag()}_{step}_{self.cfg.final_step}"
        if suffix_force:
            return f"{tag}_force"
        return tag

    def final_force_tag(self) -> str:
        return f"{self.base_tag()}_{self.cfg.final_step}_force"

    def output_tag(self) -> str:
        return f"{self.base_tag(include_double_underscore=False)}_{self.cfg.final_step}"

    def force_file(self, step: int) -> Path:
        return self.input_dir("force") / f"Force_{self.stepped_tag(step)}.txt"

    def boundary_file(self) -> Path:
        return self.input_dir("area") / f"Boundary_nodes_{self.final_force_tag()}.txt"

    def connectivity_file(self) -> Path:
        return self.input_dir("strain") / f"Lattice_connect_{self.stepped_tag(0)}.txt"
    

    def p1_connectivity_file(self) -> Path:
        p1_root = Path("lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp")
        size = "" if self.cfg.lattice_size == 64 else f"{self.cfg.lattice_size}_"
        tag = (
            f"srand_{self.cfg.seed}_{size}"
            f"{self.cfg.num_centers}_{self.cfg.num_dipoles}_1.00_"
            f"{self.cfg.tol_str}_{self.cfg.kappa_str}__{self.cfg.rest_length_str}_"
            f"{self.cfg.mu_str}_{self.cfg.mu_c_str}_0_{self.cfg.final_step}_force"
        )
        return self.cfg.base_dir / p1_root / "txt" / "strain" / f"Lattice_connect_{tag}.txt"

    def position_file(self, step: int) -> Path:
        return self.input_dir("strain") / f"Lattice_{self.stepped_tag(step)}.txt"

    def strain_file(self, step: int) -> Path:
        return self.input_dir("strain") / f"Strain_Lattice_{self.stepped_tag(step)}.txt"

    def rest_lengths_file(self, step: int) -> Path:
        return self.input_dir("rlen") / f"rlen_{self.stepped_tag(step, suffix_force=False)}.txt"

    def dipole_file(self) -> Path:
        return self.input_dir("dip_nodes") / f"Unique_dipole_nodes_{self.final_force_tag()}.txt"

    def output_file(self, subdir: str, prefix: str, ext: str = "txt") -> Path:
        return self.output_dir(subdir) / f"{prefix}_{self.output_tag()}.{ext}"

    def means_file(self, prefix: str) -> Path:
        return self.means_dir() / f"{prefix}_{self.output_tag()}.txt"

# =============================================================================
# Utility functions
# =============================================================================


def require_file(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    return path


def angle_between_points(x1: float, y1: float, x2: float, y2: float) -> float:
    return float(np.arctan2(y2 - y1, x2 - x1))


def loadtxt(path: Path) -> np.ndarray:
    print(f"Reading: {path}")
    return np.loadtxt(require_file(path))


def save_table(path: Path, columns: tuple | np.ndarray, header: str, fmt: str | tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, np.column_stack(columns), header=header, fmt=fmt)
    print(f"Wrote: {path}")


# =============================================================================
# Data loading
# =============================================================================


@dataclass
class SimulationData:
    force_all: np.ndarray
    force_magnitude: np.ndarray
    boundary_nodes: np.ndarray
    connectivity: np.ndarray
    connectivity_p1: np.ndarray | None
    positions: np.ndarray
    xpos: np.ndarray
    ypos: np.ndarray
    bond_lengths: np.ndarray
    rest_lengths: np.ndarray
    strain: np.ndarray
    dipole_nodes: np.ndarray
    dipole_centers: np.ndarray


def load_force_steps(paths: PathBuilder) -> tuple[np.ndarray, np.ndarray]:
    forces = [loadtxt(paths.force_file(step)) for step in paths.cfg.force_steps]
    force_all = -np.asarray(forces)  # Stored file is dE/dx; physical force is -dE/dx.
    force_magnitude = np.linalg.norm(force_all[:, :, :2], axis=2)
    return force_all, force_magnitude


def load_boundary_nodes(paths: PathBuilder) -> np.ndarray:
    data = loadtxt(paths.boundary_file())
    return data[:, 1].astype(int)


def load_connectivity(paths: PathBuilder) -> tuple[np.ndarray, np.ndarray | None]:
    conn = np.rint(loadtxt(paths.connectivity_file())[:, 3:9]).astype(int)
    p1_path = paths.p1_connectivity_file()
    conn_p1 = None
    if p1_path.exists():
        conn_p1 = np.rint(loadtxt(p1_path)[:, 3:9]).astype(int)
    return conn, conn_p1


def load_positions_strains_and_rest_lengths(paths: PathBuilder) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    positions = []
    bond_lengths = []
    rest_lengths = []

    for step in paths.cfg.position_steps:
        positions.append(loadtxt(paths.position_file(step)))
        bond_lengths.append(loadtxt(paths.strain_file(step))[:, 1:7])
        rest_lengths.append(loadtxt(paths.rest_lengths_file(step)))

    positions_arr = np.asarray(positions)
    bond_lengths_arr = np.asarray(bond_lengths)
    rest_lengths_arr = np.asarray(rest_lengths)
    strain = np.where(bond_lengths_arr != 0, bond_lengths_arr - rest_lengths_arr, np.nan)
    return positions_arr, bond_lengths_arr, rest_lengths_arr, strain


def load_dipole_nodes(paths: PathBuilder) -> tuple[np.ndarray, np.ndarray]:
    data = loadtxt(paths.dipole_file())
    dipole_nodes = data[:, 1].astype(int)
    dipole_centers = dipole_nodes[::7]
    return dipole_nodes, dipole_centers


def load_simulation_data(paths: PathBuilder) -> SimulationData:
    force_all, force_mag = load_force_steps(paths)
    boundary_nodes = load_boundary_nodes(paths)
    connectivity, connectivity_p1 = load_connectivity(paths)
    positions, bond_lengths, rest_lengths, strain = load_positions_strains_and_rest_lengths(paths)
    dipole_nodes, dipole_centers = load_dipole_nodes(paths)

    return SimulationData(
        force_all=force_all,
        force_magnitude=force_mag,
        boundary_nodes=boundary_nodes,
        connectivity=connectivity,
        connectivity_p1=connectivity_p1,
        positions=positions,
        xpos=positions[:, :, 1],
        ypos=positions[:, :, 2],
        bond_lengths=bond_lengths,
        rest_lengths=rest_lengths,
        strain=strain,
        dipole_nodes=dipole_nodes,
        dipole_centers=dipole_centers,
    )


# =============================================================================
# Force calculations
# =============================================================================


def compute_stretching_forces(
    positions: np.ndarray,
    connectivity: np.ndarray,
    strain: np.ndarray,
    mu_tension: float,
    mu_compression: float,
    unconnected_value: int,
) -> np.ndarray:
    """Compute node forces from bond strains for every saved step.

    Returns
    -------
    forces : np.ndarray, shape (n_steps, n_nodes, 2)
        Reconstructed stretching/compression force on every node.
    """
    n_steps, n_nodes = strain.shape[:2]
    forces = np.zeros((n_steps, n_nodes, 2), dtype=float)

    for step in range(n_steps):
        for node in range(n_nodes):
            for bond in range(6):
                neighbor = connectivity[node, bond]
                strain_val = strain[step, node, bond]

                if neighbor == unconnected_value or np.isnan(strain_val):
                    continue

                dx = positions[step, neighbor, 1] - positions[step, node, 1]
                dy = positions[step, neighbor, 2] - positions[step, node, 2]
                bond_length = np.hypot(dx, dy)
                if bond_length == 0:
                    continue

                stiffness = mu_tension if strain_val >= 0 else mu_compression
                force_scalar = stiffness * strain_val
                forces[step, node, 0] += force_scalar * dx / bond_length
                forces[step, node, 1] += force_scalar * dy / bond_length

    return forces


def compute_square_boundary_nodes(lattice_size: int) -> np.ndarray:
    nodes = []
    for node in range(lattice_size * lattice_size):
        row = node // lattice_size
        col = node % lattice_size
        if row == 0 and node not in (0, lattice_size - 1):
            nodes.append(node)
        elif col == lattice_size - 1 and node not in (lattice_size - 1, lattice_size * lattice_size - 1):
            nodes.append(node)
        elif col == 0 and node not in (0, lattice_size * (lattice_size - 1)):
            nodes.append(node)
        elif row == lattice_size - 1 and node not in (lattice_size * (lattice_size - 1), lattice_size * lattice_size - 1):
            nodes.append(node)
    return np.asarray(nodes, dtype=int)


def compute_square_boundary_normal_forces(forces: np.ndarray, lattice_size: int) -> tuple[np.ndarray, np.ndarray]:
    """Return square-boundary nodes and outward-normal force components.

    Note: the original script used -Fx on the top boundary. Here the outward-normal
    component on the top boundary is -Fy, consistent with bottom using +Fy.
    """
    boundary_nodes = []
    normal_forces = []

    for node in range(lattice_size * lattice_size):
        row = node // lattice_size
        col = node % lattice_size

        if row == 0 and node not in (0, lattice_size - 1):
            boundary_nodes.append(node)
            normal_forces.append(forces[:, node, 1])
        elif col == lattice_size - 1 and node not in (lattice_size - 1, lattice_size * lattice_size - 1):
            boundary_nodes.append(node)
            normal_forces.append(-forces[:, node, 0])
        elif col == 0 and node not in (0, lattice_size * (lattice_size - 1)):
            boundary_nodes.append(node)
            normal_forces.append(forces[:, node, 0])
        elif row == lattice_size - 1 and node not in (lattice_size * (lattice_size - 1), lattice_size * lattice_size - 1):
            boundary_nodes.append(node)
            normal_forces.append(-forces[:, node, 1])

    return np.asarray(boundary_nodes, dtype=int), np.asarray(normal_forces).T


def compute_bending_forces(
    final_positions: np.ndarray,
    connectivity: np.ndarray,
    kappa: float,
    unconnected_value: int,
) -> np.ndarray:
    """Estimate bending-force contribution on every node.

    This preserves the three-term bending-force logic from the research script,
    but removes global variables and repeated temporary initialization.
    """
    n_nodes = len(final_positions)
    forces_by_direction = np.zeros((n_nodes, 3, 2), dtype=float)
    x = final_positions[:, 1]
    y = final_positions[:, 2]

    def add_if_safe(node: int, direction: int, force_x: float, force_y: float) -> None:
        if np.isfinite(force_x) and np.isfinite(force_y):
            forces_by_direction[node, direction, 0] += force_x
            forces_by_direction[node, direction, 1] += force_y

    for node in range(n_nodes):
        for direction in range(3):
            opposite = (direction + 3) % 6

            # Case 1: node has both opposite neighbors.
            n1 = connectivity[node, direction]
            n2 = connectivity[node, opposite]
            if n1 != unconnected_value and n2 != unconnected_value:
                delxij = -(x[node] - x[n1])
                delyij = -(y[node] - y[n1])
                delxik = -(x[node] - x[n2])
                delyik = -(y[node] - y[n2])
                rij = np.hypot(delxij, delyij)
                rik = np.hypot(delxik, delyik)

                if rij > 0 and rik > 0:
                    cross = -(delxij * delyik - delyij * delxik)
                    f = abs(cross) / (rij * rik)
                    if abs(cross) > 0 and f < 1:
                        dfdx = ((delyik - delyij) * np.sign(cross) / (rij * rik) - f * delxij / rij**2 - f * delxik / rik**2)
                        dfdy = ((delxij - delxik) * np.sign(cross) / (rij * rik) - f * delyij / rij**2 - f * delyik / rik**2)
                        prefactor = kappa * f / np.sqrt(1 - f**2)
                        add_if_safe(node, direction, prefactor * dfdx, prefactor * dfdy)

            # Case 2: node-neighbor-next neighbor in same direction.
            n1 = connectivity[node, direction]
            if n1 != unconnected_value:
                n2 = connectivity[n1, direction]
                if n2 != unconnected_value:
                    delxij = -(x[node] - x[n1])
                    delyij = -(y[node] - y[n1])
                    delxmj = x[n1] - x[n2]
                    delymj = y[n1] - y[n2]
                    rij = np.hypot(delxij, delyij)
                    rmj = np.hypot(delxmj, delymj)

                    if rij > 0 and rmj > 0:
                        cross = delxij * delymj - delyij * delxmj
                        f = abs(cross) / (rij * rmj)
                        if abs(cross) > 0 and f < 1:
                            dfdx = -delymj * np.sign(cross) / (rij * rmj) - delxij * f / rij**2
                            dfdy = delxmj * np.sign(cross) / (rij * rmj) - delyij * f / rij**2
                            prefactor = kappa * f / np.sqrt(1 - f**2)
                            add_if_safe(node, direction, prefactor * dfdx, prefactor * dfdy)

            # Case 3: node-neighbor-next neighbor in opposite direction.
            n1 = connectivity[node, opposite]
            if n1 != unconnected_value:
                n2 = connectivity[n1, opposite]
                if n2 != unconnected_value:
                    delxik = -(x[node] - x[n1])
                    delyik = -(y[node] - y[n1])
                    delxlk = x[n1] - x[n2]
                    delylk = y[n1] - y[n2]
                    rik = np.hypot(delxik, delyik)
                    rlk = np.hypot(delxlk, delylk)

                    if rik > 0 and rlk > 0:
                        cross = -(delyik * delxlk - delxik * delylk)
                        f = abs(cross) / (rik * rlk)
                        if abs(cross) > 0 and f < 1:
                            dfdx = -delylk * np.sign(cross) / (rik * rlk) - delxik * f / rik**2
                            dfdy = delxlk * np.sign(cross) / (rik * rlk) - delyik * f / rik**2
                            prefactor = kappa * f / np.sqrt(1 - f**2)
                            add_if_safe(node, direction, prefactor * dfdx, prefactor * dfdy)

    return forces_by_direction.sum(axis=1)


def radial_vectors(node_positions: np.ndarray, center_positions: np.ndarray) -> np.ndarray:
    """Vectors from nodes toward the center, shape (n_steps, n_nodes, 2)."""
    return center_positions[:, None, :] - node_positions


def project_along_vectors(forces: np.ndarray, vectors: np.ndarray) -> np.ndarray:
    """Elementwise projection of force vectors along the supplied vectors."""
    norms = np.linalg.norm(vectors, axis=2)
    with np.errstate(divide="ignore", invalid="ignore"):
        projection = np.sum(forces * vectors, axis=2) / norms
    return np.nan_to_num(projection)


# =============================================================================
# Dipole moments and geometry features
# =============================================================================


def compute_boundary_dipole_moment(boundary_forces: np.ndarray, boundary_radial_vectors: np.ndarray) -> dict[str, np.ndarray | float]:
    final_forces = boundary_forces[-1]
    final_vectors = boundary_radial_vectors[-1]

    xx = final_forces[:, 0] * final_vectors[:, 0]
    yy = final_forces[:, 1] * final_vectors[:, 1]
    xy = final_forces[:, 0] * final_vectors[:, 1]
    yx = final_forces[:, 1] * final_vectors[:, 0]

    trace_per_node = xx + yy
    return {
        "xx_per_node": xx,
        "yy_per_node": yy,
        "xy_per_node": xy,
        "yx_per_node": yx,
        "trace_per_node": trace_per_node,
        "trace": float(np.sum(trace_per_node)),
        "total": float(np.sum(xx + yy + xy + yx)),
    }


def compute_local_dipole_moment(
    cfg: AnalysisConfig,
    data: SimulationData,
    bending_forces: np.ndarray,
) -> float:
    """Compute local force-dipole moment from outer dipole nodes."""
    total = 0.0

    for start in range(0, len(data.dipole_nodes), 7):
        center = data.dipole_nodes[start]
        outer_nodes = data.dipole_nodes[start + 1 : start + 7]
        dip_moment_force_mode = 0.0

        for node in outer_nodes:
            spring_vec = np.array(
                [
                    data.xpos[-1, node] - data.xpos[-1, center],
                    data.ypos[-1, node] - data.ypos[-1, center],
                ]
            )
            spring_norm = np.linalg.norm(spring_vec)
            if spring_norm == 0:
                continue

            if cfg.run_type == "force":
                angle_val = angle_between_points(
                    data.xpos[-1, center], data.ypos[-1, center], data.xpos[-1, node], data.ypos[-1, node]
                )
                imposed_force = cfg.num_force_steps * cfg.force_value * np.array([np.cos(angle_val), np.sin(angle_val)])
                dip_moment_force_mode += float(np.dot(imposed_force, spring_vec))
                continue

            node_force = np.zeros(2, dtype=float)
            angle_spring = angle_between_points(data.xpos[-1, node], data.ypos[-1, node], data.xpos[-1, center], data.ypos[-1, center])

            for bond in range(6):
                neighbor = data.connectivity[node, bond]
                strain_val = data.strain[-1, node, bond]
                if neighbor == cfg.unconnected_value or np.isnan(strain_val):
                    continue

                if neighbor == center:
                    active_strain = (cfg.rest_length + strain_val) - 1.0
                else:
                    active_strain = strain_val

                dx = data.xpos[-1, neighbor] - data.xpos[-1, node]
                dy = data.ypos[-1, neighbor] - data.ypos[-1, node]
                length = np.hypot(dx, dy)
                if length == 0:
                    continue

                stiffness = cfg.mu_tension if strain_val >= 0 else cfg.mu_compression
                force_scalar = stiffness * active_strain
                node_force += force_scalar * np.array([dx, dy]) / length

                # Kept for comparability with the old projection-based validation.
                _ = force_scalar * np.cos(angle_spring)

            node_force += bending_forces[node]
            total += float(np.dot(node_force, spring_vec))

        if cfg.run_type == "force":
            total += dip_moment_force_mode

    return total


def compute_pairwise_distances(points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if len(points) < 2:
        return np.empty((0, 2), dtype=int), np.empty(0, dtype=float)

    deltas = points[:, None, :] - points[None, :, :]
    distance_matrix = np.linalg.norm(deltas, axis=2)
    i_idx, j_idx = np.triu_indices(len(points), k=1)
    pair_indices = np.column_stack((i_idx, j_idx))
    return pair_indices, distance_matrix[i_idx, j_idx]


def compute_radius_of_gyration(points: np.ndarray) -> float:
    if len(points) == 0:
        return np.nan
    center_of_mass = points.mean(axis=0)
    return float(np.sqrt(np.mean(np.sum((points - center_of_mass) ** 2, axis=1))))


def compute_convex_hull_area(points: np.ndarray) -> float:
    if len(points) < 3:
        return 0.0
    hull = ConvexHull(points)
    return float(hull.area)


# =============================================================================
# Stress calculations
# =============================================================================


def compute_node_stress_tensor(
    node: int,
    final_positions: np.ndarray,
    connectivity: np.ndarray,
    final_strain: np.ndarray,
    bending_forces: np.ndarray,
    center_position_initial: np.ndarray,
    mu: float,
    unconnected_value: int,
) -> np.ndarray:
    """Compute 2x2 stress tensor at one node from stretching + bending terms."""
    sigma = np.zeros((2, 2), dtype=float)
    x = final_positions[:, 1]
    y = final_positions[:, 2]

    for bond in range(6):
        neighbor = connectivity[node, bond]
        strain_val = final_strain[node, bond]
        if neighbor == unconnected_value or np.isnan(strain_val):
            continue

        dx = x[neighbor] - x[node]
        dy = y[neighbor] - y[node]
        length = np.hypot(dx, dy)
        if length == 0:
            continue

        force_scalar = mu * strain_val
        force_vec = force_scalar * np.array([dx, dy]) / length
        bond_vec = np.array([dx, dy])
        sigma += np.outer(force_vec, bond_vec)

    for direction in range(3):
        n1 = connectivity[node, direction]
        n2 = connectivity[node, direction + 3]
        if n1 == unconnected_value or n2 == unconnected_value:
            continue

        rel1 = np.array([x[n1], y[n1]]) - center_position_initial
        rel2 = np.array([x[node], y[node]]) - center_position_initial
        rel3 = np.array([x[n2], y[n2]]) - center_position_initial
        f1 = bending_forces[n1]
        f2 = bending_forces[node]
        f3 = bending_forces[n2]

        # 1/3 follows common per-atom stress partitioning for angle-like interactions.
        sigma += (np.outer(f1, rel1) + np.outer(f2, rel2) + np.outer(f3, rel3)) / 3.0

    cell_area = np.sqrt(3.0) * 0.5
    return sigma / cell_area


def compute_radial_stress_from_tensor(stress_tensor: np.ndarray, node_position: np.ndarray, center_position: np.ndarray) -> float:
    radial = node_position - center_position
    norm = np.linalg.norm(radial)
    if norm == 0:
        return 0.0
    e_r = radial / norm
    return float(e_r @ stress_tensor @ e_r)


def compute_radial_stress_for_nodes(
    nodes: np.ndarray,
    data: SimulationData,
    bending_forces: np.ndarray,
    cfg: AnalysisConfig,
) -> np.ndarray:
    center_initial = np.array([data.xpos[0, cfg.center_node], data.ypos[0, cfg.center_node]])
    center_final = np.array([data.xpos[-1, cfg.center_node], data.ypos[-1, cfg.center_node]])
    final_positions = data.positions[-1]
    radial_stresses = []

    for node in nodes:
        stress = compute_node_stress_tensor(
            node=node,
            final_positions=final_positions,
            connectivity=data.connectivity,
            final_strain=data.strain[-1],
            bending_forces=bending_forces,
            center_position_initial=center_initial,
            mu=cfg.mu_tension,
            unconnected_value=cfg.unconnected_value,
        )
        node_position = np.array([data.xpos[-1, node], data.ypos[-1, node]])
        radial_stresses.append(compute_radial_stress_from_tensor(stress, node_position, center_final))

    return np.asarray(radial_stresses)


def find_nodes_in_ring(
    xpos: np.ndarray,
    ypos: np.ndarray,
    center_position: np.ndarray,
    radius: float,
    half_width: float,
) -> np.ndarray:
    distances = np.hypot(xpos - center_position[0], ypos - center_position[1])
    return np.where((distances > radius - half_width) & (distances < radius + half_width))[0]


def compute_ring_averaged_stress(
    data: SimulationData,
    bending_forces: np.ndarray,
    cfg: AnalysisConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    radii = np.arange(1, cfg.outer_radius + 1)
    mean_stress = []
    std_stress = []
    center_position = np.array([data.xpos[0, cfg.center_node], data.ypos[0, cfg.center_node]])

    for radius in radii:
        nodes = find_nodes_in_ring(data.xpos[-1], data.ypos[-1], center_position, radius, cfg.ring_width)
        if len(nodes) == 0:
            mean_stress.append(np.nan)
            std_stress.append(np.nan)
            continue
        stresses = compute_radial_stress_for_nodes(nodes, data, bending_forces, cfg)
        mean_stress.append(np.nanmean(stresses))
        std_stress.append(np.nanstd(stresses))

    return radii, np.asarray(mean_stress), np.asarray(std_stress)


# =============================================================================
# Plotting
# =============================================================================


def save_boundary_force_plots(plot_dir: Path, cfg: AnalysisConfig, radial_forces: np.ndarray) -> None:
    final_forces = radial_forces[-1]
    total = np.sum(final_forces)
    normalized = final_forces / total if total != 0 else final_forces
    tag = f"{cfg.num_centers}_{cfg.num_dipoles}_{cfg.pbond_str}_{cfg.kappa_str}_{cfg.tol_str}_{cfg.rest_length_str}_{cfg.mu_str}_{cfg.mu_c_str}"

    def scatter(values: np.ndarray, ylabel: str, filename: str, ylim: tuple[float, float] | None = None) -> None:
        plt.figure(figsize=(8, 5), dpi=300)
        plt.scatter(np.arange(len(values)), values, marker=".", s=250, alpha=0.5, label=cfg.pbond_str)
        plt.xlabel("Boundary-node index", fontsize=13)
        plt.ylabel(ylabel, fontsize=13)
        if ylim is not None:
            plt.ylim(*ylim)
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(plot_dir / filename)
        plt.close()

    scatter(final_forces, "Radial boundary force", f"bndry_force_scatter_{tag}.png")
    scatter(normalized, "Radial boundary force / total", f"bndry_force_scatter_normalized_{tag}.png")
    scatter(normalized, "Radial boundary force / total", f"bndry_force_scatter_normalized_limited_{tag}.png", ylim=(-0.1, 0.1))

    for sign_name, selector in [("positive", final_forces > 0), ("negative", final_forces < 0)]:
        values = np.abs(final_forces[selector])
        if len(values) == 0:
            continue
        plt.figure(figsize=(8, 5), dpi=300)
        plt.scatter(np.arange(len(values)), values, marker=".", s=250, alpha=0.5, label=cfg.pbond_str)
        plt.yscale("log")
        plt.xlabel("Boundary-node index", fontsize=13)
        plt.ylabel("|Radial boundary force|", fontsize=13)
        plt.title(f"{sign_name.capitalize()} radial boundary forces", fontsize=13)
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(plot_dir / f"bndry_force_scatter_{sign_name}_{tag}.png")
        plt.close()


def save_stress_force_plot(path: Path, stress: np.ndarray, force: np.ndarray, xlabel: str, ylabel: str) -> None:
    valid = np.isfinite(stress) & np.isfinite(force)
    if np.count_nonzero(valid) < 2:
        return
    slope, intercept = np.polyfit(stress[valid], force[valid], 1)
    plt.figure(figsize=(8, 5), dpi=300)
    plt.scatter(stress, force, marker="o", s=50, alpha=0.25)
    xfit = np.linspace(np.nanmin(stress), np.nanmax(stress), 200)
    plt.plot(xfit, slope * xfit + intercept, label=f"{slope:.4f} x + {intercept:.2e}")
    plt.xlabel(xlabel, fontsize=13)
    plt.ylabel(ylabel, fontsize=13)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(path)
    plt.close()


def save_ring_stress_plots(plot_dir: Path, cfg: AnalysisConfig, radii: np.ndarray, mean: np.ndarray, std: np.ndarray) -> None:
    tag = f"{cfg.lattice_size}_{cfg.num_centers}_{cfg.num_dipoles}_{cfg.pbond_str}_{cfg.tol_str}_{cfg.kappa_str}_{cfg.rest_length_str}_{cfg.mu_str}_{cfg.mu_c_str}_{cfg.final_step}"
    valid = np.isfinite(mean)
    if not np.any(valid):
        return

    xfit = np.linspace(radii[valid][0], radii[valid][-1], 200)
    yfit = mean[valid][0] / xfit**2

    plt.figure(figsize=(8, 5), dpi=300)
    plt.errorbar(radii, mean, yerr=std, linestyle="None", marker=".", ms=10, alpha=0.6, label=f"p={cfg.pbond_str}")
    plt.plot(xfit, yfit, label=r"$1/r^2$")
    plt.xlabel("Radial distance", fontsize=13)
    plt.ylabel("Mean radial stress", fontsize=13)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(plot_dir / f"radial_stress_vs_radial_dist_{tag}.png")
    plt.close()

    plt.figure(figsize=(8, 5), dpi=300)
    plt.errorbar(radii, mean, yerr=std, linestyle="None", marker=".", ms=10, alpha=0.6, label=f"p={cfg.pbond_str}")
    plt.plot(xfit, yfit, label=r"$1/r^2$")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Radial distance", fontsize=13)
    plt.ylabel("Mean radial stress", fontsize=13)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(plot_dir / f"radial_stress_vs_radial_dist_log_{tag}.png")
    plt.close()


# =============================================================================
# Pipeline
# =============================================================================


def run_analysis(cfg: AnalysisConfig) -> None:
    paths = PathBuilder(cfg)
    data = load_simulation_data(paths)

    final_positions = data.positions[-1]
    center_positions = np.column_stack((data.xpos[:, cfg.center_node], data.ypos[:, cfg.center_node]))

    # -------------------------------------------------------------------------
    # Boundary forces from stored force files
    # -------------------------------------------------------------------------
    boundary_forces = data.force_all[:, data.boundary_nodes, :2]
    
    # Match the old code:
    # force_boundary_vec[i] was projected using radial_vec[i + 1].
    boundary_positions = np.stack((data.xpos[1:, data.boundary_nodes], data.ypos[1:, data.boundary_nodes]), axis=2)
    center_positions_force_steps = center_positions[1:]
    
    boundary_radial_vecs = radial_vectors(boundary_positions, center_positions_force_steps)
    radial_boundary_force = project_along_vectors(boundary_forces, boundary_radial_vecs)

    save_table(
        paths.output_file("area", "bndry_node_force"),
        (data.boundary_nodes, boundary_forces[-1, :, 0], boundary_forces[-1, :, 1], radial_boundary_force[-1]),
        header="# node Fx Fy radial_force",
        fmt=("%12d", "%15.7e", "%15.7e", "%15.7e"),
    )

    save_table(
        paths.output_file("area", "bndry_force"),
        (cfg.rest_length_steps, np.sum(radial_boundary_force, axis=1)),
        header="# rest_length total_radial_boundary_force",
        fmt=("%15.7f", "%15.7e"),
    )

    # -------------------------------------------------------------------------
    # Reconstruct stretching forces and square-boundary forces
    # -------------------------------------------------------------------------
    stretching_forces = compute_stretching_forces(
        positions=data.positions,
        connectivity=data.connectivity,
        strain=data.strain,
        mu_tension=cfg.mu_tension,
        mu_compression=cfg.mu_compression,
        unconnected_value=cfg.unconnected_value,
    )

    for step_idx, forces_step in enumerate(stretching_forces):
        save_table(
            paths.output_file("force", f"all_nodes_forces_{step_idx}"),
            (np.arange(cfg.lattice_size * cfg.lattice_size), forces_step[:, 0], forces_step[:, 1]),
            header="# node Fx Fy",
            fmt=("%8d", "%20.7e", "%20.7e"),
        )

    square_boundary_nodes, square_boundary_normal_forces = compute_square_boundary_normal_forces(stretching_forces, cfg.lattice_size)
    for step_idx, normal_force in enumerate(square_boundary_normal_forces):
        save_table(
            paths.output_file("force", f"bndry_nodes_force_{step_idx}"),
            (square_boundary_nodes, normal_force),
            header="# node normal_force",
            fmt=("%8d", "%20.7e"),
        )

    sum_square_boundary = square_boundary_normal_forces[1:].sum(axis=1) if len(square_boundary_normal_forces) > 1 else square_boundary_normal_forces.sum(axis=1)
    rest_steps_for_square = cfg.rest_length_steps[: len(sum_square_boundary)]
    save_table(
        paths.output_file("force", "bndry_force_sum"),
        (rest_steps_for_square, sum_square_boundary),
        header="# rest_length total_square_boundary_force",
        fmt=("%15.7f", "%20.7e"),
    )

    # -------------------------------------------------------------------------
    # Boundary and local dipole moments
    # -------------------------------------------------------------------------
    boundary_moment = compute_boundary_dipole_moment(boundary_forces, boundary_radial_vecs)
    trace_per_node = boundary_moment["trace_per_node"]
    trace_moment = float(boundary_moment["trace"])

    save_table(
        paths.output_file("area", "bndry_node_moment"),
        (
            data.boundary_nodes,
            boundary_forces[-1, :, 0],
            boundary_forces[-1, :, 1],
            boundary_radial_vecs[-1, :, 0],
            boundary_radial_vecs[-1, :, 1],
            boundary_moment["xx_per_node"],
            boundary_moment["yy_per_node"],
            trace_per_node,
        ),
        header="# node Fx Fy radial_x radial_y xx_moment yy_moment trace_moment",
        fmt=("%12d", "%15.7e", "%15.7e", "%15.7e", "%15.7e", "%15.7e", "%15.7e", "%15.7e"),
    )

    save_table(
        paths.output_file("area", "bndry_node_dipole_moment"),
        (np.array([cfg.rest_length]), np.array([trace_moment])),
        header="# rest_length far_field_dipole_trace",
        fmt=("%12.3f", "%15.7e"),
    )

    bending_forces = compute_bending_forces(
        final_positions=final_positions,
        connectivity=data.connectivity,
        kappa=cfg.kappa,
        unconnected_value=cfg.unconnected_value,
    )

    save_table(
        paths.output_file("area", "bend_force_test"),
        (np.arange(cfg.lattice_size * cfg.lattice_size), bending_forces[:, 0], bending_forces[:, 1]),
        header="# node F_bend_x F_bend_y",
        fmt=("%12d", "%15.7e", "%15.7e"),
    )

    local_moment = compute_local_dipole_moment(cfg, data, bending_forces)
    ratio = trace_moment / local_moment if local_moment != 0 else np.nan

    save_table(
        paths.output_file("area", "dipole_moment_ratio"),
        (np.array([cfg.rest_length]), np.array([local_moment]), np.array([trace_moment]), np.array([ratio])),
        header="# rest_length D_local D_far D_far_over_D_local",
        fmt=("%12.3f", "%15.7e", "%15.7e", "%15.7e"),
    )

    positive_far_field = np.sum(trace_per_node[trace_per_node > 0])
    print(f"D_local = {local_moment:.7e}")
    print(f"D_far   = {trace_moment:.7e}")
    print(f"D_far / D_local = {ratio:.7e}")
    print(f"Positive far-field contributions sum to {positive_far_field:.7e}")

    # -------------------------------------------------------------------------
    # Stress on circular boundary and inner boundary
    # -------------------------------------------------------------------------
    radial_stress_boundary = compute_radial_stress_for_nodes(data.boundary_nodes, data, bending_forces, cfg)
    save_table(
        paths.output_file("displacement", "bndry_node_radial_stress"),
        (data.boundary_nodes, radial_stress_boundary, radial_boundary_force[-1]),
        header="# node sigma_rr radial_boundary_force",
        fmt=("%12d", "%15.7e", "%15.7e"),
    )

    save_stress_force_plot(
        paths.plot_dir() / "radial_stress_vs_radial_boundary_force.png",
        radial_stress_boundary,
        radial_boundary_force[-1],
        xlabel="Radial stress",
        ylabel="Radial boundary force",
    )

    center_initial = np.array([data.xpos[0, cfg.center_node], data.ypos[0, cfg.center_node]])
    inner_nodes = find_nodes_in_ring(
        data.xpos[-1], data.ypos[-1], center_initial, cfg.inner_boundary_radius, cfg.inner_boundary_width
    )

    save_table(
        paths.output_file("area", "inner_bndry_nodes"),
        (inner_nodes,),
        header="# node",
        fmt=("%12d",),
    )

    inner_forces = data.force_all[:, inner_nodes, :2]
    
    # Match the old code's force-step to position-step alignment:
    # force index i uses position/radial vector index i + 1.
    inner_positions = np.stack((data.xpos[1:, inner_nodes], data.ypos[1:, inner_nodes]), axis=2)
    
    inner_radial_vecs = radial_vectors(inner_positions, center_positions_force_steps)
    inner_radial_forces = project_along_vectors(inner_forces, inner_radial_vecs)

    inner_stress = compute_radial_stress_for_nodes(inner_nodes, data, bending_forces, cfg)

    save_stress_force_plot(
        paths.plot_dir() / "inner_radial_stress_vs_radial_boundary_force.png",
        inner_stress,
        inner_radial_forces[-1],
        xlabel="Inner radial stress",
        ylabel="Inner radial force",
    )

    save_table(
        paths.output_file("area", "total_inner_stress"),
        (np.array([len(inner_stress)]), np.array([np.sum(inner_stress)])),
        header="# n_inner_nodes total_sigma_rr",
        fmt=("%12d", "%15.7e"),
    )

    # -------------------------------------------------------------------------
    # Ring-averaged radial stress decay
    # -------------------------------------------------------------------------
    radii, mean_stress, std_stress = compute_ring_averaged_stress(data, bending_forces, cfg)
    save_table(
        paths.output_file("area", "stress_ring"),
        (radii, mean_stress, std_stress),
        header="# radius mean_sigma_rr std_sigma_rr",
        fmt=("%12d", "%15.7e", "%15.7e"),
    )
    save_ring_stress_plots(paths.plot_dir("radial_stress"), cfg, radii, mean_stress, std_stress)

    # -------------------------------------------------------------------------
    # Dipole-center geometric features
    # -------------------------------------------------------------------------
    centers = np.column_stack((data.xpos[0, data.dipole_centers], data.ypos[0, data.dipole_centers]))
    pair_indices, pairwise_distances = compute_pairwise_distances(centers)

    if cfg.pbond != 1 and cfg.num_centers != 1 and len(pairwise_distances) > 0:
        node_pairs = np.column_stack((data.dipole_centers[pair_indices[:, 0]], data.dipole_centers[pair_indices[:, 1]]))
        save_table(
            paths.means_file("center_dist"),
            (node_pairs[:, 0], node_pairs[:, 1], pairwise_distances),
            header="# center1 center2 distance",
            fmt=("%12d", "%12d", "%15.7f"),
        )

        radius_gyration = compute_radius_of_gyration(centers)
        hull_area = compute_convex_hull_area(centers)
        save_table(
            paths.means_file("centers_quant"),
            (np.array([radius_gyration]), np.array([hull_area])),
            header="# radius_gyration convex_hull_area",
            fmt=("%15.7f", "%15.7f"),
        )
        print(f"Mean center distance = {pairwise_distances.mean():.7f}")
        print(f"Std center distance  = {pairwise_distances.std():.7f}")
        print(f"Radius of gyration  = {radius_gyration:.7f}")
        print(f"Convex hull area    = {hull_area:.7f}")

    # -------------------------------------------------------------------------
    # Plots
    # -------------------------------------------------------------------------
    save_boundary_force_plots(paths.plot_dir(), cfg, radial_boundary_force)


# =============================================================================
# CLI
# =============================================================================


def parse_args() -> AnalysisConfig:
    parser = argparse.ArgumentParser(description="Analyze force, stress, and geometry features from elastic-network simulations.")
    parser.add_argument("--base-dir", type=Path, default=Path("/home/abhinav/david"), help="Base simulation-output directory.")
    parser.add_argument("--lattice-size", type=int, default=64)
    parser.add_argument("--num-centers", type=int, default=5)
    parser.add_argument("--num-force-steps", type=int, default=10)
    parser.add_argument("--pbond", type=float, default=0.55)
    parser.add_argument("--tolerance", type=float, default=1.0e-7)
    parser.add_argument("--kappa", type=float, default=1.0e-6)
    parser.add_argument("--rest-length", type=float, default=0.9)
    parser.add_argument("--mu", type=float, default=1.0, dest="mu_tension")
    parser.add_argument("--mu-compression", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=112)
    parser.add_argument("--network-id", type=str, default="111111,101111")
    parser.add_argument("--ranseed", type=str, default="667720,601210")
    parser.add_argument(
        "--run-type",
        choices=[
            "default",
            "inner_fixed",
            "radial",
            "inner_radial",
            "inner_radial_arp",
            "radial_only",
            "cluster",
            "cluster_bash",
            "cluster_hex",
            "cluster_random_hex",
            "cluster_radial",
            "force",
        ],
        default="cluster_radial",
    )
    parser.add_argument("--only-last", action="store_true")
    parser.add_argument("--only-last-num", type=int, default=2)
    parser.add_argument("--force-value", type=float, default=0.01)
    parser.add_argument("--dont-throw-disconnected-dipoles", action="store_true")
    parser.add_argument("--inner-boundary-radius", type=float, default=13.0)
    parser.add_argument("--inner-boundary-width", type=float, default=0.5)
    parser.add_argument("--outer-radius", type=int, default=25)
    parser.add_argument("--ring-width", type=float, default=0.5)
    args = parser.parse_args()

    return AnalysisConfig(
        base_dir=args.base_dir,
        lattice_size=args.lattice_size,
        num_centers=args.num_centers,
        num_force_steps=args.num_force_steps,
        pbond=args.pbond,
        tolerance=args.tolerance,
        kappa=args.kappa,
        rest_length=args.rest_length,
        mu_tension=args.mu_tension,
        mu_compression=args.mu_compression,
        seed=args.seed,
        network_id=args.network_id,
        ranseed=args.ranseed,
        run_type=args.run_type,
        only_last=args.only_last,
        only_last_num=args.only_last_num,
        force_value=args.force_value,
        dont_throw_disconnected_dipoles=args.dont_throw_disconnected_dipoles,
        inner_boundary_radius=args.inner_boundary_radius,
        inner_boundary_width=args.inner_boundary_width,
        outer_radius=args.outer_radius,
        ring_width=args.ring_width,
    )


def main() -> None:
    cfg = parse_args()
    run_analysis(cfg)


if __name__ == "__main__":
    main()
