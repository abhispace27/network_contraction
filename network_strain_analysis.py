#!/usr/bin/env python3
"""Analyze strain, stress, and displacement in circular elastic-network simulations.

This script is a cleaned-up, GitHub-ready version of a research analysis script. It
keeps the original file-naming convention, but wraps the workflow in functions,
uses pathlib/argparse, and keeps the vectorized NumPy operations explicit.

Typical usage
-------------
python network_strain_analysis.py \
    --base-dir /home/abhinav/david \
    --lattice-size 64 \
    --num-centers 5 \
    --pbond 0.55 \
    --kappa 1e-5 \
    --srand 112 \
    --diff-network 111111,101111 \
    --plot-boundary-forces

Notes
-----
The script assumes the same folder and filename structure produced by the
simulation code. It does not change the physics/math of the analysis; the main
changes are organization, safer path handling, command-line options, and removal
of inactive notebook-style blocks.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib import cm
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import griddata
from scipy.optimize import curve_fit


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class AnalysisConfig:
    base_dir: Path
    lattice_size: int = 64
    num_centers: int = 5
    num_steps: int = 11
    pbond: float = 0.55
    mu: float = 1.0
    mu_c: float = 1.0
    tol: float = 1e-7
    kappa: float = 1e-5
    rlen: float = 0.9
    force_val: float = 0.01
    srand: int | None = 112
    ranseed: str = "667720,601210"
    diff_network: str | None = "111111,101111"
    pbd: bool = False
    buckle: bool = False
    anisotropic: bool = False
    force_mode: bool = False
    folder_mode: str = "cluster_radial"
    only_last: bool = False
    only_last_num: int = 2
    plot_strain: bool = False
    plot_bending_heatmap: bool = False
    plot_boundary_forces: bool = True
    color_strain: bool = False
    movie: bool = False
    compute_displacement: bool = True

    @property
    def num_points(self) -> int:
        return self.lattice_size * self.lattice_size

    @property
    def num_dipoles(self) -> int:
        return self.num_centers if self.anisotropic else self.num_centers * 6

    @property
    def strain_threshold(self) -> float:
        return 1e-3 if self.pbond == 1 else 1e-6

    @property
    def pbond_txt(self) -> str:
        return f"{self.pbond:.2f}"

    @property
    def tol_txt(self) -> str:
        return f"{self.tol:.2e}"

    @property
    def kappa_txt(self) -> str:
        return "0.00e+00" if self.kappa == 0 else f"{self.kappa:.2e}"

    @property
    def rlen_txt(self) -> str:
        return f"{self.rlen:.4f}"

    @property
    def force_txt(self) -> str:
        # Original script used extra precision in force mode.
        return f"{self.force_val:.9f}" if self.force_mode else f"{self.force_val:.4f}"

    @property
    def mu_txt(self) -> str:
        return f"{self.mu:.4f}"

    @property
    def mu_c_txt(self) -> str:
        return f"{self.mu_c:.4f}"

    @property
    def value_txt(self) -> str:
        return self.force_txt if self.force_mode else self.rlen_txt

    @property
    def inner_radius(self) -> float:
        return 12.0

    @property
    def outer_radius(self) -> float:
        return 50.0 if self.lattice_size == 128 else 25.0

    @property
    def center_node(self) -> int:
        return 8383 if self.lattice_size == 128 else 2080

    @property
    def unconnected_value(self) -> int:
        return -1 if self.lattice_size == 128 else 9999

    @property
    def force_steps(self) -> list[int]:
        if not self.only_last:
            return list(range(self.num_steps))
        return list(range(self.only_last_num))


@dataclass
class LoadedData:
    xpos: np.ndarray
    ypos: np.ndarray
    strain_raw: np.ndarray
    rlen_raw: np.ndarray
    conn_data: np.ndarray
    conn_node: np.ndarray
    boundary_nodes_circle: np.ndarray
    outer_nodes: np.ndarray
    inner_nodes: np.ndarray
    dipole_nodes: np.ndarray
    bending_energy_node: np.ndarray


@dataclass
class BondPlotData:
    x_all: list[np.ndarray]
    y_all: list[np.ndarray]
    x_plot: list[np.ndarray]
    y_plot: list[np.ndarray]
    square_boundary_nodes: np.ndarray
    node_bond_plot: np.ndarray
    flat_plot_index: np.ndarray


@dataclass
class StrainData:
    strain_all: np.ndarray
    strain_circle: np.ndarray
    circle_nodes: np.ndarray


# -----------------------------------------------------------------------------
# Path helpers
# -----------------------------------------------------------------------------


def kappa_folder(kappa: float) -> str:
    mapping = {
        0.0: "kappa2_0",
        2e-7: "kappa2_2e-7",
        5e-7: "kappa2_5e-7",
        1e-6: "kappa2_e-6",
        2e-6: "kappa2_2e-6",
        5e-6: "kappa2_5e-6",
        1e-5: "kappa2_e-5",
        2e-5: "kappa2_2e-5",
        5e-5: "kappa2_5e-5",
        1e-4: "kappa2_e-4",
        2e-4: "kappa2_2e-4",
        5e-4: "kappa2_5e-4",
        1e-3: "kappa2_e-3",
        2e-3: "kappa2_2e-3",
        5e-3: "kappa2_5e-3",
        1e-2: "kappa2_e-2",
    }
    try:
        return mapping[float(kappa)]
    except KeyError as exc:
        allowed = ", ".join(str(k) for k in sorted(mapping))
        raise ValueError(f"Unsupported kappa={kappa}. Allowed values: {allowed}") from exc


def folder_root(cfg: AnalysisConfig) -> str:
    folder_map = {
        "radial": "lattice_nopbd_bndry_all_clamp_restlength_circle_radial",
        "inner_radial": "lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial",
        "inner_radial_arp": "lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp",
        "radial_only": "lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_only_radial_arp",
        "inner_radial_arp_bash": "lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash",
        "cluster": "lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_cluster",
        "cluster_new": "cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash",
        "cluster_radial": "cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_radial",
        "hex": "cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_hex",
        "hex_rand": "cluster_download/lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_bash_random_dipole_hex_bonds",
        "force": "lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_arp_force",
        "anisotropic": "lattice_nopbd_bndry_all_clamp_restlength_circle_arp_anisotropic",
        "default_inner_fixed": "lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed",
        "default": "lattice_nopbd_bndry_all_clamp_restlength_circle",
    }

    if cfg.force_mode:
        return folder_map["force"]
    if cfg.anisotropic:
        return folder_map["anisotropic"]
    if cfg.folder_mode in folder_map:
        return folder_map[cfg.folder_mode]
    raise ValueError(f"Unknown folder mode: {cfg.folder_mode}")


def simulation_folder(cfg: AnalysisConfig) -> Path:
    root = folder_root(cfg)

    if cfg.pbond == 1:
        if cfg.lattice_size == 64 and cfg.srand == 300:
            root = "lattice_nopbd_bndry_all_clamp_restlength_circle_inner_fixed_radial_large"
        return cfg.base_dir / root

    parts: list[str] = [root, kappa_folder(cfg.kappa), cfg.ranseed]
    if cfg.srand is not None:
        parts.append(str(cfg.srand))
    if cfg.diff_network:
        parts.append(cfg.diff_network)
    return cfg.base_dir.joinpath(*parts)


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def loadtxt_checked(path: Path, *, dtype=float) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    return np.loadtxt(path, dtype=dtype)


def file_stem(cfg: AnalysisConfig, step: int, final_step: int | None = None) -> str:
    final = cfg.num_steps - 1 if final_step is None else final_step
    return (
        f"{cfg.num_centers}_{cfg.num_dipoles}_{cfg.pbond_txt}_{cfg.tol_txt}_"
        f"{cfg.kappa_txt}__{cfg.value_txt}_{cfg.mu_txt}_{cfg.mu_c_txt}_{step}_{final}"
    )


def prefixed_name(cfg: AnalysisConfig, prefix: str, step: int, suffix: str = "") -> str:
    stem = file_stem(cfg, step)
    if cfg.srand is not None and (cfg.pbond == 1 or cfg.force_mode):
        stem = f"srand_{cfg.srand}_{stem}"
    if cfg.lattice_size != 64:
        stem = f"{cfg.lattice_size}_{stem}" if cfg.pbond < 1 else f"srand_{cfg.srand}_{cfg.lattice_size}_{stem}"
    return f"{prefix}{stem}{suffix}"


def area_name(cfg: AnalysisConfig, label: str) -> str:
    suffix = "_force.txt"
    final = cfg.num_steps - 1
    value = cfg.value_txt
    stem = (
        f"{cfg.num_centers}_{cfg.num_dipoles}_{cfg.pbond_txt}_{cfg.tol_txt}_"
        f"{cfg.kappa_txt}__{value}_{cfg.mu_txt}_{cfg.mu_c_txt}_{final}{suffix}"
    )
    if cfg.lattice_size != 64:
        stem = f"{cfg.lattice_size}_{stem}" if cfg.pbond < 1 else f"srand_{cfg.srand}_{cfg.lattice_size}_{stem}"
    elif cfg.srand is not None and (cfg.pbond == 1 or cfg.force_mode):
        stem = f"srand_{cfg.srand}_{stem}"
    return f"{label}_{stem}"


# -----------------------------------------------------------------------------
# Geometry/vectorization helpers
# -----------------------------------------------------------------------------


def angle(x1: np.ndarray, y1: np.ndarray, x2: np.ndarray, y2: np.ndarray) -> np.ndarray:
    return np.arctan2(y2 - y1, x2 - x1)


def flatten_segments(x0: np.ndarray, y0: np.ndarray, x1: np.ndarray, y1: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return np.column_stack((x0, x1)).ravel(), np.column_stack((y0, y1)).ravel()


def line_segments_from_flat_xy(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.stack(
        (np.column_stack((x[0::2], y[0::2])), np.column_stack((x[1::2], y[1::2]))),
        axis=1,
    )


def smooth_circular_vectors(values: np.ndarray) -> np.ndarray:
    return (np.roll(values, 1) + values + np.roll(values, -1)) / 3.0


def build_bond_plot_arrays(
    xpos: np.ndarray,
    ypos: np.ndarray,
    conn_node: np.ndarray,
    outer_nodes: np.ndarray,
    cfg: AnalysisConfig,
) -> BondPlotData:
    n_steps, num_points = xpos.shape
    nodes = np.arange(num_points)
    rows = nodes // cfg.lattice_size
    cols = nodes % cfg.lattice_size
    boundary_nodes = nodes[(rows == 0) | (rows == cfg.lattice_size - 1) | (cols == 0) | (cols == cfg.lattice_size - 1)]

    valid = conn_node[:, :3] != cfg.unconnected_value
    src, bond = np.where(valid)
    dst = conn_node[src, bond]

    plot_mask = ~(np.isin(src, outer_nodes) | np.isin(dst, outer_nodes))
    x_all: list[np.ndarray] = []
    y_all: list[np.ndarray] = []
    x_plot: list[np.ndarray] = []
    y_plot: list[np.ndarray] = []

    for step in range(n_steps):
        x0 = xpos[step, src]
        y0 = ypos[step, src]
        x1 = xpos[step, dst].copy()
        y1 = ypos[step, dst]

        if cfg.pbd:
            left = cols[src] == 0
            right = cols[src] == cfg.lattice_size - 1
            x1 = np.where(left & ((x1 - x0) > cfg.lattice_size / 2), x1 - cfg.lattice_size, x1)
            x1 = np.where(right & ((x0 - x1) > cfg.lattice_size / 2), x1 + cfg.lattice_size, x1)

        x_flat, y_flat = flatten_segments(x0, y0, x1, y1)
        xp_flat, yp_flat = flatten_segments(x0[plot_mask], y0[plot_mask], x1[plot_mask], y1[plot_mask])
        x_all.append(x_flat)
        y_all.append(y_flat)
        x_plot.append(xp_flat)
        y_plot.append(yp_flat)

    node_bond_single = np.column_stack((src[plot_mask], bond[plot_mask])).astype(int)
    node_bond_plot = np.broadcast_to(node_bond_single, (n_steps, *node_bond_single.shape)).copy()
    flat_plot_index = node_bond_single[:, 0] * 3 + node_bond_single[:, 1]

    return BondPlotData(x_all, y_all, x_plot, y_plot, boundary_nodes, node_bond_plot, flat_plot_index)


# -----------------------------------------------------------------------------
# Loading
# -----------------------------------------------------------------------------


def data_prefixes(cfg: AnalysisConfig, folder: Path) -> dict[str, Path]:
    strain_dir = folder / "txt" / "strain"
    rlen_dir = folder / "txt" / "rlen"
    energy_dir = folder / "energy" / "node"

    lattice_prefix = "Lattice_"
    strain_prefix = "Strain_Lattice_"
    rlen_prefix = "rlen_"
    energy_prefix = "Lattice_node_"

    if cfg.srand is not None and cfg.pbond == 1:
        lattice_prefix += f"srand_{cfg.srand}_"
        strain_prefix += f"srand_{cfg.srand}_"
        rlen_prefix += f"srand_{cfg.srand}_"
        energy_prefix += f"srand_{cfg.srand}_"

    if cfg.lattice_size != 64 and not (cfg.pbond == 1 and cfg.srand is not None):
        lattice_prefix += f"{cfg.lattice_size}_"
        strain_prefix += f"{cfg.lattice_size}_"
        rlen_prefix += f"{cfg.lattice_size}_"

    return {
        "pos": strain_dir / lattice_prefix,
        "strain": strain_dir / strain_prefix,
        "rlen": rlen_dir / rlen_prefix,
        "energy": energy_dir / energy_prefix,
    }


def load_position_strain_rlen(cfg: AnalysisConfig, folder: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    prefixes = data_prefixes(cfg, folder)
    positions: list[np.ndarray] = []
    strain_rows: list[np.ndarray] = []
    rlen_rows: list[np.ndarray] = []

    for i in cfg.force_steps:
        step_for_filename = i
        if cfg.only_last and i == cfg.only_last_num - 1 and cfg.srand not in {123, 124}:
            # Fixed from the original condition: (srand != 123 or srand != 124) was always true.
            step_for_filename = cfg.num_steps - 1

        stem = file_stem(cfg, step_for_filename)
        pos_path = Path(str(prefixes["pos"]) + f"{stem}_force.txt")
        strain_path = Path(str(prefixes["strain"]) + f"{stem}_force.txt")
        rlen_path = Path(str(prefixes["rlen"]) + f"{stem}.txt")

        print(f"Loading position file: {pos_path}")
        positions.append(loadtxt_checked(pos_path))
        strain_rows.append(loadtxt_checked(strain_path))
        rlen_rows.append(loadtxt_checked(rlen_path))

    return np.asarray(positions), np.asarray(strain_rows), np.asarray(rlen_rows)


def connection_path(cfg: AnalysisConfig, folder: Path) -> Path:
    prefix = "Lattice_connect_"
    if cfg.lattice_size != 64 and cfg.pbond != 1:
        prefix += f"{cfg.lattice_size}_"
    elif cfg.lattice_size != 64 and cfg.pbond == 1:
        prefix += f"srand_{cfg.srand}_{cfg.lattice_size}_"
    elif cfg.srand is not None and cfg.pbond == 1:
        prefix += f"srand_{cfg.srand}_"
    elif cfg.force_mode:
        prefix += f"srand_{cfg.srand}_"

    return folder / "txt" / "strain" / f"{prefix}{file_stem(cfg, 0)}_force.txt"


def load_area_nodes(cfg: AnalysisConfig, folder: Path, label: str) -> np.ndarray:
    path = folder / "txt" / "area" / area_name(cfg, label)
    print(f"Loading {label.lower()} file: {path}")
    data = loadtxt_checked(path)
    return data[:, 1].astype(int) if data.ndim == 2 and data.shape[1] > 1 else data.astype(int)


def load_dipole_nodes(cfg: AnalysisConfig, folder: Path) -> np.ndarray:
    path = folder / "txt" / "dip_nodes" / area_name(cfg, "Unique_dipole_nodes")
    print(f"Loading dipole-node file: {path}")
    data = loadtxt_checked(path)
    if cfg.num_centers <= 0:
        return np.array([0], dtype=int)
    return np.rint(data[:, 1]).astype(int)


def load_bending_energy(cfg: AnalysisConfig, folder: Path) -> np.ndarray:
    prefixes = data_prefixes(cfg, folder)
    step = cfg.only_last_num - 1 if cfg.only_last else cfg.num_steps - 1
    stem = (
        f"{cfg.num_centers}_{cfg.num_dipoles}_{cfg.pbond_txt}_{cfg.tol_txt}_{cfg.kappa_txt}_"
        f"{cfg.rlen_txt}_{cfg.mu_txt}_{cfg.mu_c_txt}_{step}_{cfg.num_steps - 1}_force.txt"
    )
    path = Path(str(prefixes["energy"]) + stem)
    print(f"Loading bending-energy file: {path}")
    data = loadtxt_checked(path)
    return np.sum(data[:, 7:], axis=1)


def load_all_data(cfg: AnalysisConfig, folder: Path) -> LoadedData:
    all_data, strain_raw, rlen_raw = load_position_strain_rlen(cfg, folder)
    xpos = all_data[:, :, 1]
    ypos = all_data[:, :, 2]

    conn_data = loadtxt_checked(connection_path(cfg, folder))
    conn_node = np.rint(conn_data[:, 3:9]).astype(int)

    boundary_nodes = load_area_nodes(cfg, folder, "Boundary_nodes")
    outer_nodes = load_area_nodes(cfg, folder, "Outer_nodes")
    inner_nodes = load_area_nodes(cfg, folder, "Inner_nodes") if cfg.diff_network else np.array([], dtype=int)
    dipole_nodes = load_dipole_nodes(cfg, folder)
    bending_energy = load_bending_energy(cfg, folder)

    # Original analysis excludes the circular boundary from the outer-node mask.
    outer_nodes = np.setdiff1d(outer_nodes, boundary_nodes, assume_unique=False)

    return LoadedData(
        xpos=xpos,
        ypos=ypos,
        strain_raw=strain_raw,
        rlen_raw=rlen_raw[:, :, :3],
        conn_data=conn_data,
        conn_node=conn_node,
        boundary_nodes_circle=boundary_nodes,
        outer_nodes=outer_nodes,
        inner_nodes=inner_nodes,
        dipole_nodes=dipole_nodes,
        bending_energy_node=bending_energy,
    )


# -----------------------------------------------------------------------------
# Analysis
# -----------------------------------------------------------------------------


def compute_strain(data: LoadedData, plot_data: BondPlotData, cfg: AnalysisConfig) -> StrainData:
    n_steps = len(cfg.force_steps)
    rlen_flat = data.rlen_raw.reshape(n_steps, 3 * cfg.num_points)
    bond_len = data.strain_raw[:, :, 1:4].reshape(n_steps, 3 * cfg.num_points)

    all_nodes = np.arange(cfg.num_points)
    circle_nodes = np.setdiff1d(all_nodes, data.outer_nodes, assume_unique=False)

    strain_circle = bond_len[:, plot_data.flat_plot_index] - rlen_flat[:, plot_data.flat_plot_index]

    strain_all = np.zeros_like(bond_len)
    valid = bond_len != 0
    strain_all[valid] = (bond_len - rlen_flat)[valid]

    # Preserve the original convention that the first saved state is the reference.
    strain_all[0, :] = 0.0
    strain_circle[0, :] = 0.0

    if cfg.buckle:
        strain_all = np.where(strain_all < 0, strain_all * cfg.mu_c, strain_all)
        strain_circle = np.where(strain_circle < 0, strain_circle * cfg.mu_c, strain_circle)

    return StrainData(strain_all=strain_all, strain_circle=strain_circle, circle_nodes=circle_nodes)


def compute_inner_boundary_stress(
    data: LoadedData,
    plot_data: BondPlotData,
    strain_data: StrainData,
    cfg: AnalysisConfig,
) -> dict[str, list[np.ndarray] | np.ndarray]:
    node_cross: list[np.ndarray] = []
    bond_cross: list[np.ndarray] = []
    angle_cross_list: list[np.ndarray] = []
    strain_x_list: list[np.ndarray] = []
    strain_y_list: list[np.ndarray] = []
    strain_cross_list: list[np.ndarray] = []
    mid_x_list: list[np.ndarray] = []
    mid_y_list: list[np.ndarray] = []
    radial_force_list: list[np.ndarray] = []

    strain_circle = strain_data.strain_circle.copy()

    for step in range(plot_data.node_bond_plot.shape[0]):
        nodes = plot_data.node_bond_plot[step, :, 0]
        bonds = plot_data.node_bond_plot[step, :, 1]
        connected = data.conn_node[nodes, bonds]

        dr_node = np.hypot(data.xpos[step, nodes] - data.xpos[step, cfg.center_node], data.ypos[step, nodes] - data.ypos[step, cfg.center_node])
        dr_conn = np.hypot(data.xpos[step, connected] - data.xpos[step, cfg.center_node], data.ypos[step, connected] - data.ypos[step, cfg.center_node])

        cross_mask = ((dr_node < cfg.inner_radius) & (dr_conn > cfg.inner_radius)) | ((dr_node > cfg.inner_radius) & (dr_conn < cfg.inner_radius))
        cross_idx = np.flatnonzero(cross_mask)

        # Original sign convention: stretched bonds give negative stress.
        strain_circle[step, cross_idx] *= -1
        strain_cross = strain_circle[step, cross_idx]

        nodes_cross = nodes[cross_idx]
        bonds_cross = bonds[cross_idx]
        conn_cross = connected[cross_idx]

        theta = angle(data.xpos[step, nodes_cross], data.ypos[step, nodes_cross], data.xpos[step, conn_cross], data.ypos[step, conn_cross])
        strain_x = strain_cross * np.cos(theta)
        strain_y = strain_cross * np.sin(theta)

        mid_x = 0.5 * (data.xpos[step, nodes_cross] + data.xpos[step, conn_cross])
        mid_y = 0.5 * (data.ypos[step, nodes_cross] + data.ypos[step, conn_cross])
        radial_norm = np.hypot(mid_x, mid_y)
        radial_force = np.divide(strain_x * mid_x + strain_y * mid_y, radial_norm, out=np.zeros_like(strain_cross), where=radial_norm != 0)

        node_cross.append(nodes_cross)
        bond_cross.append(bonds_cross)
        angle_cross_list.append(theta)
        strain_x_list.append(strain_x)
        strain_y_list.append(strain_y)
        strain_cross_list.append(strain_cross)
        mid_x_list.append(mid_x)
        mid_y_list.append(mid_y)
        radial_force_list.append(radial_force)

    total_radial_stress = np.array([np.sum(values) / (2 * np.pi * cfg.inner_radius) for values in radial_force_list])

    return {
        "node_cross": node_cross,
        "bond_cross": bond_cross,
        "angle_cross": angle_cross_list,
        "strain_x": strain_x_list,
        "strain_y": strain_y_list,
        "strain": strain_cross_list,
        "mid_x": mid_x_list,
        "mid_y": mid_y_list,
        "radial_force": radial_force_list,
        "total_radial_stress": total_radial_stress,
    }


def write_stress_outputs(stress: dict[str, list[np.ndarray] | np.ndarray], cfg: AnalysisConfig, folder: Path) -> None:
    out_dir = folder / "txt" / "force"
    out_dir.mkdir(parents=True, exist_ok=True)

    last = -1
    base = (
        f"{cfg.srand}_{cfg.num_centers}_{cfg.num_dipoles}_{cfg.pbond_txt}_{cfg.tol_txt}_"
        f"{cfg.kappa_txt}_{cfg.rlen_txt}_{cfg.mu_txt}_{cfg.mu_c_txt}_{cfg.num_steps - 1}"
    )
    if cfg.lattice_size != 64:
        base = f"{cfg.srand}_{cfg.lattice_size}_{cfg.num_centers}_{cfg.num_dipoles}_{cfg.pbond_txt}_{cfg.tol_txt}_{cfg.kappa_txt}_{cfg.rlen_txt}_{cfg.mu_txt}_{cfg.mu_c_txt}_{cfg.num_steps - 1}"

    individual_path = out_dir / f"inner_bndry_bond_radial_force_{base}.txt"
    total_path = out_dir / f"inner_bndry_bond_radial_force_total_{base}.txt"

    table = np.column_stack(
        (
            stress["node_cross"][last],
            stress["bond_cross"][last],
            stress["angle_cross"][last],
            stress["mid_x"][last],
            stress["mid_y"][last],
            stress["strain_x"][last],
            stress["strain_y"][last],
            stress["strain"][last],
            stress["radial_force"][last],
        )
    )
    np.savetxt(
        individual_path,
        table,
        header="node bond bond_angle pos_x pos_y strain_x strain_y strain radial_strain",
        fmt=("%12d", "%12d", "%15.7e", "%15.7e", "%15.7e", "%15.7e", "%15.7e", "%15.7e", "%15.7e"),
    )

    total = np.column_stack((len(stress["node_cross"][last]), stress["total_radial_stress"][last]))
    np.savetxt(total_path, total, header="total_nodes radial_stress", fmt=("%12d", "%15.7e"))
    print(f"Wrote stress outputs:\n  {individual_path}\n  {total_path}")


def write_coordination_output(data: LoadedData, cfg: AnalysisConfig, folder: Path) -> None:
    all_nodes = np.arange(cfg.num_points)
    non_dip_nodes = np.setdiff1d(all_nodes, data.dipole_nodes)
    non_dip_nodes = np.setdiff1d(non_dip_nodes, data.outer_nodes)

    coord_num = data.conn_data[:, -1]
    mean_coord = np.mean(coord_num[non_dip_nodes])
    std_coord = np.std(coord_num[non_dip_nodes])

    out_dir = folder / "txt" / "area"
    out_dir.mkdir(parents=True, exist_ok=True)
    name = f"mean_coord_num_{cfg.srand}_{cfg.num_centers}_{cfg.num_dipoles}_{cfg.pbond_txt}_{cfg.tol_txt}_{cfg.kappa_txt}_{cfg.rlen_txt}_{cfg.mu_txt}_{cfg.mu_c_txt}_{cfg.num_steps - 1}.txt"
    if cfg.lattice_size != 64:
        name = f"mean_coord_num_{cfg.srand}_{cfg.lattice_size}_{cfg.num_centers}_{cfg.num_dipoles}_{cfg.pbond_txt}_{cfg.tol_txt}_{cfg.kappa_txt}_{cfg.rlen_txt}_{cfg.mu_txt}_{cfg.mu_c_txt}_{cfg.num_steps - 1}.txt"

    out_path = out_dir / name
    np.savetxt(out_path, np.column_stack((len(non_dip_nodes), mean_coord, std_coord)), header="total_nodes mean_z std_z", fmt=("%12d", "%15.5e", "%15.5e"))
    print(f"Wrote coordination output: {out_path}")


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------


def sci_notation(value: float, pos: int | None = None) -> str:
    if value == 0:
        return "0"
    exponent = int(np.floor(np.log10(abs(value))))
    mantissa = value / 10**exponent
    return rf"${mantissa:.1f}\times10^{{{exponent}}}$"


def strain_colormap() -> ListedColormap:
    red = cm.get_cmap("Reds_r", 512)
    blue = cm.get_cmap("Blues", 512)
    red_small = ListedColormap(red(np.linspace(0.0, 0.7, 128)))
    blue_small = ListedColormap(blue(np.linspace(0.3, 1.0, 128)))
    colors = np.vstack((red_small(np.linspace(0, 1, 128)), blue_small(np.linspace(0, 1, 128))))
    return ListedColormap(colors, name="red-blue-strain")


def transformed_strain(values: np.ndarray, cfg: AnalysisConfig, log_scale: bool = False) -> np.ndarray:
    out = values.copy()
    if log_scale:
        pos = out > 0
        neg = out < 0
        out[pos] = -np.log(out[pos])
        out[neg] = np.log(np.abs(out[neg]))
    return np.where(out < 0, out * cfg.mu_c, out)


def node_style(data: LoadedData, cfg: AnalysisConfig, highlight_ring: Iterable[int] | None = None) -> tuple[np.ndarray, np.ndarray]:
    colors = np.full(cfg.num_points, "grey", dtype=object)
    sizes = np.full(cfg.num_points, 0.2)

    if data.inner_nodes.size:
        sizes[data.inner_nodes] = 3.0 if cfg.lattice_size == 128 else 10.0
        colors[data.inner_nodes] = "k"

    if highlight_ring is not None:
        ring = np.asarray(list(highlight_ring), dtype=int)
        sizes[ring] = 20.0
        colors[ring] = "m"

    sizes[data.dipole_nodes] = 200.0 if cfg.movie else (10.0 if cfg.lattice_size == 128 else 50.0)
    colors[data.dipole_nodes] = "g"
    sizes[data.outer_nodes] = 0.0
    colors[data.outer_nodes] = "w"
    sizes[data.boundary_nodes_circle] = 20.0
    colors[data.boundary_nodes_circle] = "m"
    return sizes, colors


def make_strain_collection(step: int, plot_data: BondPlotData, strain_data: StrainData, cfg: AnalysisConfig) -> tuple[LineCollection, float]:
    raw = strain_data.strain_circle[step]
    threshold = abs(cfg.strain_threshold)
    segments = line_segments_from_flat_xy(plot_data.x_plot[step], plot_data.y_plot[step])

    widths = np.full(raw.shape, 0.5)
    stretch = raw > threshold
    compress = raw < -threshold
    neutral = np.abs(raw) < threshold
    widths[stretch | compress] *= 4

    if cfg.color_strain:
        values = transformed_strain(raw, cfg)
        vlim = np.nanmax(np.abs(values)) if values.size else 1.0
        collection = LineCollection(segments, lw=widths, array=values, cmap=strain_colormap(), norm=plt.Normalize(-vlim, vlim))
    else:
        colors = np.full(raw.shape, "k", dtype=object)
        colors[stretch] = "b"
        colors[compress] = "r"
        colors[neutral] = "k"
        collection = LineCollection(segments, lw=widths, colors=colors)

    return collection, threshold


def format_network_axis(ax: plt.Axes, cfg: AnalysisConfig) -> None:
    if cfg.movie:
        ax.set_xlim([30, 50])
        ax.set_ylim([35, 50])
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.tick_params(axis="both", which="both", left=False, bottom=False)


def add_colorbar(fig: plt.Figure, ax: plt.Axes, collection: LineCollection) -> None:
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(collection, ax=ax, cax=cax, format=mticker.FuncFormatter(sci_notation))
    cbar.ax.tick_params(labelsize=18)
    cbar.locator = mticker.MaxNLocator(nbins=3)
    cbar.update_ticks()


def plot_name(cfg: AnalysisConfig, folder: Path, prefix: str, step: int, threshold: float, extra: str = "") -> Path:
    plot_dir = folder / "png" / "strain"
    plot_dir.mkdir(parents=True, exist_ok=True)
    if cfg.movie:
        return plot_dir / f"img{step:03d}.jpg"

    stem = (
        f"srand_{cfg.srand}_{cfg.lattice_size}_{cfg.num_centers}_{cfg.num_dipoles}_{cfg.pbond_txt}_"
        f"{cfg.mu_txt}_{cfg.mu_c_txt}_{cfg.kappa_txt}_{cfg.tol_txt}_{cfg.value_txt}_{threshold:.2e}_{cfg.lattice_size}_{step}_{cfg.num_steps - 1}"
    )
    return plot_dir / f"{prefix}_{stem}{extra}.png"


def plot_circle_strain(step: int, data: LoadedData, plot_data: BondPlotData, strain_data: StrainData, cfg: AnalysisConfig, folder: Path) -> None:
    collection, threshold = make_strain_collection(step, plot_data, strain_data, cfg)
    sizes, colors = node_style(data, cfg)

    fig, ax = plt.subplots(figsize=(15, 9), dpi=100)
    ax.add_collection(collection)
    ax.scatter(data.xpos[step], data.ypos[step], s=sizes, marker="o", color=colors, alpha=0.5)
    if cfg.color_strain:
        add_colorbar(fig, ax, collection)
    format_network_axis(ax, cfg)
    fig.subplots_adjust(right=0.85, left=0.05, top=0.95, bottom=0.05)

    prefix = "circle_colored_strain" if cfg.color_strain else "circle_binary_strain"
    output = plot_name(cfg, folder, prefix, step, threshold)
    fig.savefig(output, dpi=100 if cfg.movie else 600)
    plt.close(fig)
    print(f"Saved strain plot: {output}")


def plot_bending_heatmap(step: int, data: LoadedData, plot_data: BondPlotData, strain_data: StrainData, cfg: AnalysisConfig, folder: Path) -> None:
    collection, threshold = make_strain_collection(step, plot_data, strain_data, cfg)

    xi = np.linspace(np.min(data.xpos[-1]), np.max(data.xpos[-1]), 300)
    yi = np.linspace(np.min(data.ypos[-1]), np.max(data.ypos[-1]), 300)
    xi, yi = np.meshgrid(xi, yi)
    zi = griddata((data.xpos[-1], data.ypos[-1]), data.bending_energy_node, (xi, yi), method="cubic")

    fig, ax = plt.subplots(figsize=(15, 9), dpi=300)
    cmap_orange = LinearSegmentedColormap.from_list("fade_white_orange", [(1, 1, 1), (1, 0.5, 0)])
    heatmap = ax.contourf(xi, yi, zi, levels=100, cmap=cmap_orange, vmin=0)
    cbar = fig.colorbar(heatmap, ax=ax)
    cbar.set_label("Bending Energy", fontsize=18)
    cbar.ax.yaxis.set_major_formatter(mticker.FuncFormatter(sci_notation))

    ax.add_collection(collection)
    ring_rad = 25 if cfg.lattice_size == 128 else 12
    dist = np.hypot(data.xpos[0] - data.xpos[0, cfg.center_node], data.ypos[0] - data.ypos[0, cfg.center_node])
    ring_nodes = np.flatnonzero((dist > ring_rad - 0.5) & (dist < ring_rad + 0.5))
    sizes, colors = node_style(data, cfg, highlight_ring=ring_nodes)
    ax.scatter(data.xpos[step], data.ypos[step], s=sizes, marker="o", color=colors, alpha=0.5)

    format_network_axis(ax, cfg)
    output = plot_name(cfg, folder, "circle_heatmap", step, threshold)
    fig.savefig(output, dpi=300)
    plt.close(fig)
    print(f"Saved bending heatmap: {output}")


def boundary_force_path(cfg: AnalysisConfig, folder: Path) -> Path:
    stem = (
        f"{cfg.num_centers}_{cfg.num_dipoles}_{cfg.pbond_txt}_{cfg.tol_txt}_{cfg.kappa_txt}_"
        f"{cfg.rlen_txt}_{cfg.mu_txt}_{cfg.mu_c_txt}_{cfg.num_steps - 1}"
    )
    if cfg.force_mode or (cfg.srand is not None and cfg.pbond == 1 and cfg.lattice_size == 64):
        stem = f"srand_{cfg.srand}_{stem}"
    elif cfg.srand is not None and cfg.lattice_size != 64:
        stem = f"srand_{cfg.srand}_{cfg.lattice_size}_{stem}"
    return folder / "txt" / "area" / f"bndry_node_force_{stem}.txt"


def boundary_quiver_data(force_data: np.ndarray, data: LoadedData, step: int, cfg: AnalysisConfig) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    x_boundary = data.xpos[step, data.boundary_nodes_circle]
    y_boundary = data.ypos[step, data.boundary_nodes_circle]
    theta = np.mod(np.arctan2(y_boundary - np.mean(y_boundary), x_boundary - np.mean(x_boundary)), 2 * np.pi)
    sort_idx = np.argsort(theta)
    sample_idx = np.arange(0, len(sort_idx), 2)

    fx = force_data[:, 1]
    fy = force_data[:, 2]
    fx_plot = smooth_circular_vectors(fx[sort_idx])[sample_idx]
    fy_plot = smooth_circular_vectors(fy[sort_idx])[sample_idx]

    if cfg.num_centers in {5, 15} and cfg.pbond == 0.55:
        quiver_scale = 1e-5
    elif cfg.num_centers == 20 and cfg.pbond == 0.55:
        quiver_scale = 1e-6
    else:
        quiver_scale = 100 * np.mean(np.hypot(fx_plot, fy_plot))

    if cfg.num_centers == 5 and cfg.kappa == 1e-6:
        max_mag = 1e-6
    elif cfg.num_centers == 5 and cfg.kappa == 1e-5:
        max_mag = 2e-6
    elif cfg.num_centers == 20:
        max_mag = 1.5e-7
    else:
        max_mag = 2e-6

    mag = np.hypot(fx_plot, fy_plot)
    clip = np.minimum(1.0, max_mag / np.maximum(mag, 1e-30))
    return x_boundary[sort_idx][sample_idx], y_boundary[sort_idx][sample_idx], fx_plot * clip, fy_plot * clip, quiver_scale


def plot_boundary_forces(step: int, data: LoadedData, plot_data: BondPlotData, strain_data: StrainData, cfg: AnalysisConfig, folder: Path) -> None:
    force_path = boundary_force_path(cfg, folder)
    print(f"Loading boundary-force file: {force_path}")
    force_data = loadtxt_checked(force_path)

    collection, threshold = make_strain_collection(step, plot_data, strain_data, cfg)
    sizes, colors = node_style(data, cfg)
    xq, yq, fxq, fyq, quiver_scale = boundary_quiver_data(force_data, data, step, cfg)

    fig, ax = plt.subplots(figsize=(15, 9), dpi=100)
    ax.add_collection(collection)
    ax.scatter(data.xpos[step], data.ypos[step], s=sizes, marker="o", color=colors, alpha=0.5)
    ax.scatter(data.xpos[step, data.dipole_nodes], data.ypos[step, data.dipole_nodes], s=20, marker="o", color="g", zorder=50)
    ax.quiver(xq, yq, fxq, fyq, scale=quiver_scale, zorder=10)
    if cfg.color_strain:
        add_colorbar(fig, ax, collection)
    ax.set_axis_off()
    format_network_axis(ax, cfg)
    fig.subplots_adjust(right=0.85, left=0.05, top=0.95, bottom=0.05)

    output = plot_name(cfg, folder, "force_bndry", step, threshold, extra=f"_{quiver_scale:.4f}")
    fig.savefig(output, dpi=100 if cfg.movie else 600)
    plt.close(fig)
    print(f"Saved boundary-force plot: {output}")


# -----------------------------------------------------------------------------
# Displacement and fitting
# -----------------------------------------------------------------------------


def continuum_radial_displacement(x: np.ndarray, sigma: float, r1: float, r2: float) -> np.ndarray:
    return sigma * r1**2 / (2 * ((np.sqrt(3) / 2) * r1**2 + (np.sqrt(3) / 4) * r2**2)) * ((r2**2 - x**2) / x)


def compute_displacement_outputs(data: LoadedData, cfg: AnalysisConfig, folder: Path) -> None:
    x0 = data.xpos[0]
    y0 = data.ypos[0]
    x1 = data.xpos[-1]
    y1 = data.ypos[-1]

    init_r = np.hypot(x0 - x0[cfg.center_node], y0 - y0[cfg.center_node])
    final_r = np.hypot(x1 - x1[cfg.center_node], y1 - y1[cfg.center_node])
    radial_disp = final_r - init_r

    xlim = cfg.outer_radius + 0.5
    bins = np.arange(1, xlim, 1)
    dd = 0.5
    mask = (init_r[None, :] >= bins[:, None] - dd) & (init_r[None, :] < bins[:, None] + dd)
    counts = mask.sum(axis=1)
    disp_masked = np.where(mask, radial_disp[None, :], np.nan)
    mean_disp = np.divide(np.nansum(disp_masked, axis=1), counts, out=np.full(len(bins), np.nan), where=counts > 0)
    median_disp = np.nanmedian(disp_masked, axis=1)
    std_disp = np.nanstd(disp_masked, axis=1)

    out_dir = folder / "txt" / "displacement"
    out_dir.mkdir(parents=True, exist_ok=True)

    suffix = f"{cfg.srand}_{cfg.num_centers}_{cfg.num_dipoles}_{cfg.pbond_txt}_{cfg.tol_txt}_{cfg.kappa_txt}_{cfg.rlen_txt}_{cfg.mu_txt}_{cfg.mu_c_txt}_{cfg.num_steps - 1}.txt"
    if cfg.lattice_size != 64:
        suffix = f"srand_{cfg.srand}_{cfg.lattice_size}_{cfg.num_centers}_{cfg.num_dipoles}_{cfg.pbond_txt}_{cfg.tol_txt}_{cfg.kappa_txt}_{cfg.rlen_txt}_{cfg.mu_txt}_{cfg.mu_c_txt}_{cfg.num_steps - 1}.txt"

    np.savetxt(
        out_dir / f"bndry_node_radial_disp_{suffix}",
        np.column_stack((np.arange(cfg.num_points), init_r, final_r, radial_disp)),
        header="node initial_radial_dist final_radial_dist radial_disp",
        fmt=("%12d", "%15.7e", "%15.7e", "%15.7e"),
    )
    np.savetxt(
        out_dir / f"bndry_node_radial_disp_mean_{suffix}",
        np.column_stack((bins, mean_disp, median_disp, std_disp)),
        header="bin_mid mean_radial_disp median_radial_disp std_radial_disp",
        fmt=("%12d", "%15.7e", "%15.7e", "%15.7e"),
    )

    r1 = cfg.inner_radius + 1
    r2 = 48 if cfg.lattice_size != 64 else 25
    start = int(r1 - 1)
    xdata = bins[start:]
    ydata = mean_disp[start:]
    valid = np.isfinite(xdata) & np.isfinite(ydata)
    popt, _ = curve_fit(lambda x, sigma: continuum_radial_displacement(x, sigma, r1, r2), xdata[valid], ydata[valid], maxfev=10000)

    np.savetxt(
        out_dir / f"fit_radial_disp_{suffix}",
        np.column_stack((cfg.pbond, cfg.kappa, cfg.num_centers, popt[0])),
        header="pbond kappa N Sigma1_over_alpha_m",
        fmt=("%12.4f", "%10.2e", "%12d", "%15.7e"),
    )
    print(f"Fit parameter Sigma1/alpha_m: {popt[0]:.6e}")

    plot_displacement(init_r, radial_disp, bins, mean_disp, popt[0], r1, r2, cfg, folder)


def plot_displacement(init_r: np.ndarray, radial_disp: np.ndarray, bins: np.ndarray, mean_disp: np.ndarray, sigma_fit: float, r1: float, r2: float, cfg: AnalysisConfig, folder: Path) -> None:
    plot_dir = folder / "png"
    plot_dir.mkdir(parents=True, exist_ok=True)

    suffix = f"srand_{cfg.srand}_{cfg.lattice_size}_{cfg.num_centers}_{cfg.num_dipoles}_{cfg.pbond_txt}_{cfg.mu_txt}_{cfg.mu_c_txt}_{cfg.kappa_txt}_{cfg.tol_txt}_{cfg.rlen_txt}_{cfg.num_steps - 1}.png"

    fig, ax = plt.subplots(figsize=(8, 4), dpi=300)
    ax.scatter(init_r, radial_disp, alpha=0.3)
    ax.set_xlabel("Distance from center", fontsize=15)
    ax.set_ylabel("Radial displacement", fontsize=15)
    ax.set_xlim([1, cfg.outer_radius + 0.5])
    fig.subplots_adjust(top=0.92, left=0.15, bottom=0.15, right=0.98)
    fig.savefig(plot_dir / f"displacement_scatter_{suffix}")
    plt.close(fig)

    fit_x = np.linspace(r1, r2, 100)
    fig, ax = plt.subplots(figsize=(8, 4), dpi=300)
    ax.scatter(bins, mean_disp)
    ax.plot(fit_x, continuum_radial_displacement(fit_x, sigma_fit, r1, r2), label=rf"$\Sigma_1={sigma_fit:.2e}$")
    ax.set_xlabel("Distance from center", fontsize=15)
    ax.set_ylabel("Mean radial displacement", fontsize=15)
    ax.set_xlim([cfg.inner_radius + 0.5, cfg.outer_radius + 0.5])
    ax.legend()
    fig.subplots_adjust(top=0.92, left=0.17, bottom=0.15, right=0.98)
    fig.savefig(plot_dir / f"displacement_scatter_outer_region_mean_{suffix}")
    plt.close(fig)


# -----------------------------------------------------------------------------
# CLI / main
# -----------------------------------------------------------------------------


def parse_args() -> AnalysisConfig:
    parser = argparse.ArgumentParser(description="Analyze circular elastic-network simulation outputs.")
    parser.add_argument("--base-dir", type=Path, default=Path("/home/abhinav/david"))
    parser.add_argument("--lattice-size", type=int, default=64)
    parser.add_argument("--num-centers", type=int, default=5)
    parser.add_argument("--num-steps", type=int, default=11)
    parser.add_argument("--pbond", type=float, default=0.55)
    parser.add_argument("--mu", type=float, default=1.0)
    parser.add_argument("--mu-c", type=float, default=1.0)
    parser.add_argument("--tol", type=float, default=1e-7)
    parser.add_argument("--kappa", type=float, default=1e-5)
    parser.add_argument("--rlen", type=float, default=0.9)
    parser.add_argument("--force-val", type=float, default=0.01)
    parser.add_argument("--srand", type=int, default=112)
    parser.add_argument("--ranseed", default="667720,601210")
    parser.add_argument("--diff-network", default="111111,101111")
    parser.add_argument("--folder-mode", default="cluster_radial")
    parser.add_argument("--pbd", action="store_true")
    parser.add_argument("--buckle", action="store_true")
    parser.add_argument("--anisotropic", action="store_true")
    parser.add_argument("--force-mode", action="store_true")
    parser.add_argument("--only-last", action="store_true")
    parser.add_argument("--only-last-num", type=int, default=2)
    parser.add_argument("--plot-strain", action="store_true")
    parser.add_argument("--plot-bending-heatmap", action="store_true")
    parser.add_argument("--plot-boundary-forces", action="store_true")
    parser.add_argument("--no-boundary-forces", dest="plot_boundary_forces", action="store_false")
    parser.set_defaults(plot_boundary_forces=True)
    parser.add_argument("--color-strain", action="store_true")
    parser.add_argument("--movie", action="store_true")
    parser.add_argument("--no-displacement", dest="compute_displacement", action="store_false")
    parser.set_defaults(compute_displacement=True)
    args = parser.parse_args()

    mu_c = 0.1 if args.buckle else args.mu_c
    force_val = 0.01115 if args.force_mode and args.force_val == 0.01 else args.force_val

    return AnalysisConfig(
        base_dir=args.base_dir,
        lattice_size=args.lattice_size,
        num_centers=args.num_centers,
        num_steps=args.num_steps,
        pbond=args.pbond,
        mu=args.mu,
        mu_c=mu_c,
        tol=args.tol,
        kappa=args.kappa,
        rlen=args.rlen,
        force_val=force_val,
        srand=args.srand,
        ranseed=args.ranseed,
        diff_network=args.diff_network if args.diff_network.lower() != "none" else None,
        pbd=args.pbd,
        buckle=args.buckle,
        anisotropic=args.anisotropic,
        force_mode=args.force_mode,
        folder_mode=args.folder_mode,
        only_last=args.only_last,
        only_last_num=args.only_last_num,
        plot_strain=args.plot_strain,
        plot_bending_heatmap=args.plot_bending_heatmap,
        plot_boundary_forces=args.plot_boundary_forces,
        color_strain=args.color_strain,
        movie=args.movie,
        compute_displacement=args.compute_displacement,
    )


def main() -> None:
    cfg = parse_args()
    folder = simulation_folder(cfg)
    print(f"Simulation folder: {folder}")

    data = load_all_data(cfg, folder)
    plot_data = build_bond_plot_arrays(data.xpos, data.ypos, data.conn_node, data.outer_nodes, cfg)
    strain_data = compute_strain(data, plot_data, cfg)

    stress = compute_inner_boundary_stress(data, plot_data, strain_data, cfg)
    write_stress_outputs(stress, cfg, folder)
    write_coordination_output(data, cfg, folder)

    last_step = len(cfg.force_steps) - 1
    if cfg.plot_strain:
        for step in range(len(cfg.force_steps)):
            plot_circle_strain(step, data, plot_data, strain_data, cfg, folder)

    if cfg.plot_bending_heatmap:
        plot_bending_heatmap(last_step, data, plot_data, strain_data, cfg, folder)

    if cfg.plot_boundary_forces:
        plot_boundary_forces(last_step, data, plot_data, strain_data, cfg, folder)

    if cfg.compute_displacement:
        compute_displacement_outputs(data, cfg, folder)

    centers = data.dipole_nodes[::7]
    print(f"Dipole center nodes: {centers}")
    print("Analysis complete.")


if __name__ == "__main__":
    main()
