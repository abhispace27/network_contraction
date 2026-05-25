# network_contraction

## Publication

This simulation and analysis pipeline supported research that led to a peer-reviewed publication in *Soft Matter*. The code was used to study how localized contractile force dipoles transmit stress, strain, and displacement through disordered elastic networks.

If using this repository for academic purposes, please cite the associated publication:
> Kumar A, Quint DA, Dasbiswas K, Cooperative effect of local active stresses on the macroscopic contractility of elastic fiber networks, Soft Matter (2026) DOI: https://doi.org/10.1039/D5SM00803D  

# Details
Simulation and analysis pipeline for studying how localized active force dipoles generate long-range mechanical response in disordered elastic networks.

This repository contains C simulation code, SLURM/HPC automation scripts, and Python analysis/visualization pipelines for elastic-network contraction simulations. The project models a two-dimensional disordered spring network with circular boundaries, introduces contractile dipoles by changing local bond rest lengths, relaxes the network using energy minimization, and post-processes the output to quantify strain, force transmission, stress propagation, displacement fields, and effective dipole moments.

The goal is to connect microscopic network-level mechanics with coarse-grained observables such as radial stress, strain localization, boundary force transmission, and continuum-scale displacement decay.

---

## Project Overview

The workflow has four main parts:

1. **C simulation code**
   - Builds a 2D elastic lattice network.
   - Randomly removes bonds to generate disordered network architectures.
   - Places active contractile dipoles inside the network.
   - Changes local bond rest lengths to model active contraction.
   - Computes stretching, compression, bending, and angular energy contributions.
   - Relaxes the network using conjugate-gradient energy minimization.
   - Writes node positions, bond strains, forces, energies, boundary nodes, dipole nodes, and network connectivity to structured output files.

2. **SLURM batch script**
   - Automates large parameter sweeps on an HPC cluster.
   - Runs simulations over multiple random seeds, network realizations, bending stiffnesses, bond probabilities, and dipole numbers.
   - Uses SLURM array jobs with a configurable maximum number of concurrent tasks.
   - Dynamically computes the number of required simulations from the parameter arrays.
   - Compiles and runs one simulation per array task.

3. **`force_stress_analysis_refactored.py`**
   - Reconstructs force and stress observables from simulation output.
   - Computes stretching forces, bending-force contributions, radial boundary forces, force-dipole moments, and ring-averaged radial stress.
   - Extracts geometric features of dipole organization, including pairwise distances, radius of gyration, and convex-hull area.
   - Saves tabular outputs and summary plots.

4. **`network_strain_analysis.py`**
   - Analyzes strain, displacement, stress, and visualization outputs.
   - Computes bond-level strain from deformed and rest-length configurations.
   - Identifies bonds crossing an inner circular boundary.
   - Estimates radial stress from bond-level force projections.
   - Computes radial displacement profiles relative to the network center.
   - Performs radial binning and fits displacement profiles to continuum-theory predictions.
   - Generates high-resolution network visualizations using Matplotlib and `LineCollection`.

---

## Scientific Motivation

Active biological and soft-material networks often generate macroscopic deformation from localized microscopic force-generating units. Examples include cytoskeletal contraction, motor-driven fiber networks, active gels, and mechanically responsive disordered materials.

This project studies a simplified computational model of that problem:

> How do local contractile force dipoles transmit stress, strain, and displacement through a disordered elastic network?

The simulations allow controlled variation of network architecture, bending stiffness, bond occupation probability, and dipole organization. The analysis pipeline then measures how those microscopic parameters affect large-scale mechanical response.

---

## Repository Structure

```text
network_contraction/
├── Movies-Working_2D_Lamel_network_force_c...c
├── RAN.c
├── README.md
├── force_stress_analysis.py
├── network_strain_analysis.py
├── nrutil_jen.c
├── nrutil_jen.h
└── sim_script.sub
