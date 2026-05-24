# network_contraction
C code with energy minimization and python codes for analysis and visualizations

# --------------------------------------------------------------------------------------------------
# Analysis code: network_strain_analysis.y

This project contains a Python analysis pipeline for post-processing elastic-network simulation data. The code analyzes deformation patterns produced by active force dipoles in a disordered lattice network with circular boundaries.

The pipeline computes bond-level strain, identifies bonds crossing an inner circular boundary, estimates radial stress, measures node-level radial displacement, performs radial binning, fits displacement profiles to continuum-theory predictions, and generates high-resolution scientific visualizations of network deformation.

The goal of this project is to connect particle/network-level simulation outputs with coarse-grained mechanical observables such as strain localization, radial stress transmission, displacement decay, and effective contractile response.

## Key Features

- Loads and processes large simulation output files from elastic-network simulations
- Computes bond strain from deformed and rest-length configurations
- Identifies network bonds crossing a chosen circular boundary
- Estimates radial stress from bond-level force projections
- Computes radial displacement profiles relative to the network center
- Performs radial binning, averaging, and curve fitting using NumPy/SciPy
- Generates publication-quality network plots using Matplotlib and `LineCollection`
- Supports command-line configuration through `argparse`
- Uses vectorized NumPy operations where possible for efficient analysis

## Technical Skills Demonstrated

- Scientific Python programming
- NumPy-based vectorized data processing
- Simulation data analysis
- Network/lattice geometry analysis
- Scientific visualization with Matplotlib
- Curve fitting with SciPy
- Command-line workflow design
- Modular code organization using dataclasses and functions

## Example Use Case

This analysis was developed for studying how localized active stresses propagate through elastic fiber networks. The code can be used to quantify how microscopic force dipoles generate long-range strain, radial stress, and displacement fields in disordered mechanical networks.
