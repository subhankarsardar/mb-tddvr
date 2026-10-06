# MB-TDDVR: Multi-Basis Time-Dependent Discrete Variable Representation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Language: Fortran](https://img.shields.io/badge/Language-Fortran%202008+-blue.svg)](https://fortran-lang.org/)
[![OpenMP](https://img.shields.io/badge/Parallel-OpenMP-green.svg)](https://www.openmp.org/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23191203.svg)](https://doi.org/10.5281/zenodo.23191203)

MB-TDDVR is an open-source Fortran framework for high-dimensional
quantum dynamics simulations of vibronically coupled molecular
systems. It implements a modular, adaptive, and parallel version of
the Time-Dependent Discrete Variable Representation (TDDVR) method,
with coordinate-dependent DVR basis selection, multiple propagation
schemes, adaptive time-stepping, and OpenMP parallelization.

The code has been validated on two systems:

1. The Hénon–Heiles model (4D–32D) — benchmark against MCTDH and ML-MCTDH
2. The Barrelene radical cation (C₈H₈⁺/Bl⁺) — 42 vibrational modes, 6 coupled electronic states

---

## Table of Contents

- [Features](#features)
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Repository Structure](#repository-structure)
- [Input and Output Files](#input-and-output-files)
- [Basis Selection](#basis-selection)
- [Propagation Schemes](#propagation-schemes)
- [Analysis Modules](#analysis-modules)
- [Known Limitations](#known-limitations)
- [Testing](#testing)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

---

## Features

MB-TDDVR provides the following capabilities:

- **Coordinate-dependent DVR basis selection.** Hermite, Legendre,
  and sine DVR basis functions can be assigned independently to each
  vibrational mode according to its physical character.

- **Modular propagation suite.** Five propagators are implemented
  behind a common interface:
  - 2nd-order split-operator (SO-2)
  - 4th-order Suzuki–Trotter split-operator (SO-4)
  - 2nd-order Magnus expansion
  - 4th-order Runge–Kutta (RK4)
  - Adaptive split-operator with local error control

- **Adaptive time-stepping.** Automatic step-size adjustment based
  on local truncation error estimates. Reduces the number of
  propagation steps by 30–50% without loss of accuracy.

- **OpenMP parallelization.** Shared-memory parallelization of the
  mode-independent kinetic-energy operations.

- **Analysis modules.** Energy conservation monitoring, adiabatic
  population analysis, wavepacket visualization, and 1D/2D potential
  cuts along selected modes.

- **Modern Fortran.** Written in Fortran 2008+ using derived types,
  allocatable arrays, modules, and procedure interfaces.

---

## System Requirements

### Minimum Requirements

| Component | Requirement |
|---|---|
| **Operating System** | Linux (tested on Ubuntu 20.04/22.04, CentOS 7) |
| **Compiler** | gfortran ≥ 9.0 or ifort ≥ 19.0 |
| **Parallelization** | OpenMP (bundled with modern compilers) |
| **BLAS/LAPACK** | Reference or optimized (MKL, OpenBLAS) |
| **Memory** | 4 GB minimum; 8 GB recommended for Bl⁺ |
| **Disk** | 500 MB for source, build, and example output |

### Tested Configurations

The code has been tested with:

- gfortran 11.4.0 (Ubuntu 22.04)
- ifort 2021.5.0 (Intel oneAPI)

### Hardware Used for Benchmarks

- 32-core Intel Xeon Gold 6248R @ 3.00 GHz
- 128 GB RAM

---

## Installation

### Step 1: Clone the repository

```bash
git clone https://github.com/your-username/mb-tddvr.git
cd mb-tddvr



# Build with the provided script:
chmod +x compile-script.sh
./compile-script.sh

# Quick Start:
Running the Barrelene (Bl⁺) example, The default configuration in SO-TDDVR.f90 runs the 42-mode, 6-state Barrelene radical cation with the adaptive split-operator propagator.

./compile-script.sh
./tddvr_so > barrelene_run.log &

# Controlling the number of OpenMP threads
export OMP_NUM_THREADS=8
./tddvr_so

===========================================================
Input and Output Files
Main Output Files
File	Description
barrelene_probabilities.dat	Diabatic state populations vs. time
barrelene_autocorrelation.dat	Wavepacket autocorrelation function
barrelene_classical_trajectory.dat	Classical coordinates and momenta
barrelene_timing.dat	Performance and timing statistics
Module Output Directories
Directory	Contents
./energy_analysis/	Total energy, kinetic, potential, components, conservation check
./adiabatic_analysis/	Adiabatic populations, energies, transformation matrices
./wavepacket_analysis/	State populations, mode probabilities, 2D slices
./potential_cuts/	1D potential cuts along selected modes
./potential_2d_cuts/	2D potential surfaces for selected mode pairs


====================================================================================
Citation
If you use MB-TDDVR in your research, please cite:

Manuscript:

Sardar, S. MB-TDDVR: A Modular, Adaptive, and Parallel Framework
for Efficient High-Dimensional Quantum Dynamics. 2026.
Manuscript in preparation.

==========================================================
Contact
Author: Dr. Subhankar Sardar

Affiliation: Department of Chemistry, Bhatter College,
Dantan, Paschim Medinipur, West Bengal, India

Email: subhankar.iacs@gmail.com

Issues: Please report bugs or request features via the
GitHub issue tracker.
============================================================

