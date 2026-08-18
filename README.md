# Characterizing the Ion Conductive State of the α7 nAChR 

[![DOI](https://img.shields.io/badge/DOI-10.1021%2Facs.jpcb.5c08465-blue.svg)](https://pubs.acs.org/doi/10.1021/acs.jpcb.5c08465)

**Authors:** Nauman Sultan, Gisela D. Cymes, Ada Chen, Bernard Brooks, Claudio Grosman, Ana Damjanovic  
**Affiliations:** Johns Hopkins University | National Institutes of Health (NIH) | University of Illinois at Urbana-Champaign

## Overview
This repository contains the system preparation files, simulation inputs, analysis scripts, and plotting code used in the publication: **"Characterizing the Ion-Conductive State of the α7-Nicotinic Acetylcholine Receptor via Single-Channel Measurements and Molecular Dynamics Simulations"** (*The Journal of Physical Chemistry B*). 

The computational pipeline outlines the workflow, ranging from system setup and molecular dynamics simulations to kinetic modeling and data visualization. 

## Repository Structure
```text
a7-nAChR-MD-Conduction/
├── 1_system_setup/               # Starting structures and forcefield parameters
│   ├── apo/
│   ├── ligand_parameters/
│   └── modulator_bound/
├── 2_simulation_inputs/          # MD equilibration and production configs
│   ├── restrained_run/
│   └── unrestrained_run/
├── 3_analysis_scripts/           # Trajectory analysis and dimensionality reduction
├── 4_kinetic_modeling/           # Mathematical modeling of ion conduction
└── 5_figure_generation/          # MATLAB routines for manuscript figures
```

## System Setup (`1_system_setup`)
Contains the starting configurations and force field parameters required to build the simulation systems.

* **`apo/` & `modulator_bound/`**: Contains the base structural files (`.pdb`, `.psf`) for various receptor states (PDB IDs: 7EKT, 7KOX, 8V80, 8V82, 9LH5).
* **`ligand_parameters/`**: Contains the CHARMM/AMBER topologies and parameters (`.mol2`, `.prm`, `.rtf`, `.str`, `.frcmod`) for all specific ligands evaluated in the study (CLR, I33, I34, EPJ, YLR, NCT).

## Simulation Inputs (`2_simulation_inputs`)
Contains the configuration files to execute the Molecular Dynamics pipeline using AMBER and OpenMM.

* **`restrained_run/`**: Includes AMBER production scripts (`amber_production_restrained.in`, `amber_production_restrained_Ca.in`) and position restraints (`Restraint.rest`, `Restraint_Ca.rest`) used during the initial equilibration phase.
* **`unrestrained_run/`**: Includes the AMBER (`amber_Ca_restraints.in`, `amber_production_unrestrained.in`) with position restraints (`CAL_fix.rest`), and OpenMM input files  (`openmm_Ca_restraints.py`) for generating the final unrestrained production trajectories.

## Analysis Scripts (`3_analysis_scripts`)
The parsing utilities and analytical scripts responsible for extracting biophysical properties from the raw MD trajectories. 

* **`calc_electrostatics.sh`**: A shell script to compute the electrostatic profiles of the receptor.
* **`calc_hydration_num.py`**: A Python script to calculate water coordination and hydration numbers within the channel.
* **`calc_ion_density.py`**: Extracts the 3D spatial density of conducting ions.
* **`calc_ion_flux.py`**: Quantifies the rate and direction of ion permeation across the channel.
* **`calc_pca_tsne.py`**: Performs Principal Component Analysis (PCA) and t-SNE dimensionality reduction to classify receptor conformational states.
* **`calc_pore_radius.py`**: Computes the physical dimensions and bottleneck radius of the conductive pore over time.
* **`calc_rmsd_box.py`**: Tracks the Root Mean Square Deviation (RMSD) and simulation box dimensions to verify system stability.

## Kinetic Modeling (`4_kinetic_modeling`)
* **`DiffSolution.nb`**: A Wolfram Mathematica notebook containing the analytical diffusion solutions for the conduction model.
* **`double_poisson_fit.py`**: A Python script utilizing a double Poisson distribution to model the statistical kinetics of ion crossing events.

## Figure Generation (`5_figure_generation`)
MATLAB plotting routines used to generate the final visualizations for the manuscript.

* **`plot_channel_radius_comparison.m`**: Generates a comparative plot of the pore radius across different receptor states.
* **`plot_comprehensive_profiles.m`**: Aggregates multi-variable data into a comprehensive structural and energetic profile.
* **`plot_M2_residues.m`**: Visualizes the specific orientations and behaviors of the critical M2 pore-lining residues.
* **`plot_position_den.m`**: Plots the 1D/2D positional densities of the ions within the channel.

## Usage
1. Clone the repository to your local machine or cluster.
2. Ensure you have the required dependencies installed (AMBER, OpenMM, Python 3, MATLAB).
3. Execute the workflow sequentially from directories 1 through 5. 

*(Note: Raw MD trajectory files are not hosted directly on this repository. Execution of the scripts in step 3 requires the base trajectory data, which is available upon reasonable request to the corresponding authors).*

*(Note: The `results/` directory is not included in this repository. You will need to manually create an empty `results/` folder at the root of the repository to store the outputs from step 3. The MATLAB scripts in step 5 are hardcoded to read analysis data from this directory).*

## Citation
If you use these scripts, parameters, or methods in your research, please cite our paper:
```bibtex
@article{Sultan_2026,
  author = {Sultan, Nauman and Cymes, Gisela D. and Chen, Ada and Brooks, Bernard and Grosman, Claudio and Damjanovic, Ana},
  title = {Characterizing the Ion-Conductive State of the α7-Nicotinic Acetylcholine Receptor via Single-Channel Measurements and Molecular Dynamics Simulations},
  journal = {The Journal of Physical Chemistry B},
  year = {2026},
  doi = {10.1021/acs.jpcb.5c08465}
}
```

## License
This project is licensed under the MIT License - see the `LICENSE` file for details.