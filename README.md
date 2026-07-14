# pySIMAID - Simulator for Adsorption-Induced Deformation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![LAMMPS](https://img.shields.io/badge/LAMMPS-required-red.svg)](https://www.lammps.org/)

A Python-controlled workflow for simulating gas adsorption and structural deformation in metal-organic frameworks (MOFs) using a thermodynamically rigorous hybrid MD/GCMC approach.

## Overview

**pySIMAID** (Python Simulator for Adsorption-Induced Deformation) is designed to capture the coupled phenomena of gas adsorption and framework flexibility in nanoporous materials. The method combines:

- **GCMC** for efficient sampling of particle insertions/deletions at fixed chemical potential
- **NVT relaxation** for adsorbate relaxation and velocity assignment
- **NPT MD** for realistic structural dynamics with volume fluctuations
- **Metropolis acceptance** based on rigorous free energy comparisons

This approach enables accurate prediction of adsorption isotherm and framework deformation in flexible MOFs.

## Key Features

-  **Thermodynamically consistent** acceptance criterion
-  **Automatic restart** capability from checkpoints
-  **Flexible configuration** via command-line arguments
-  **Detailed statistics** tracking and logging
-  **Modular LAMMPS inputs** for easy force field customization
-  **MPI-parallel** execution support

## Scientific Motivation

Traditional GCMC simulations assume rigid frameworks, while purely MD-based approaches trquire unfeasibly large supercells. pySIMAID bridges this gap by:

1. Using GCMC to efficiently sample adsorption/desorption
2. Allowing the framework to relax and deform via MD
3. Maintaining detailed balance through Metropolis acceptance

This captures adsorption-induced deformation phenomena such as:
- Framework breathing and swelling
- Gate-opening transitions
- Cooperative adsorption effects
- Pressure-dependent structural changes

## Requirements

### Software Dependencies
- **Python 3.7+** (standard library only)
- **LAMMPS** (2020 or later) compiled with:
  - `MC` package (for GCMC)
  - `KSPACE` package (for long-range electrostatics)
  - `MOLECULE` package (for molecular systems)
  - MPI support
- **MPI** implementation (OpenMPI, MPICH, or Intel MPI)

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/pySIMAID.git
cd pySIMAID

# Verify LAMMPS installation
mpirun -np 1 lmp_mpi -help

# Make main script executable
chmod +x pysimaid.py

# Test installation
./pysimaid.py --help
```

## Quick Start

```bash
# Run a simulation at 87.3 K and 1 millibar
./pysimaid.py -T 87.3 -P 0.00101325 -n 1000

# Monitor progress
tail -f T_87.3/P_0.00101325/acceptance_stats.txt

# Check acceptance rate
grep -c "1$" T_87.3/P_0.00101325/acceptance_stats.txt
```

## File Structure

```
pySIMAID/
├── pysimaid.py                    # Main control script
├── lammps_inputs/
│   ├── equilibrate_empty.in       # Prepare + equilibrate empty framework
│   ├── gcmc_step.in              # GCMC equilibration
│   ├── nvt_step.in               # NVT relaxation (produces reference state)
│   ├── npt_step.in               # NPT MD trajectory
│   └── paircoeffs.in             # Force field parameters
├── structures/
│   └── *.data                   # Initial MOF structure
├── README.md
├── LICENSE
└── requirements.txt             
```

## Usage

### Basic Command

```bash
./pysimaid.py [OPTIONS]
```

### Command-Line Options

#### LAMMPS Execution
```bash
--lammps-exec PATH       # Path to LAMMPS executable (default: lmp_mpi)
--nprocs N              # Number of MPI processes (default: 16)
```

#### Thermodynamic Conditions
```bash
-T, --temperature TEMP   # Temperature in Kelvin (default: 87.3)
-P, --pressure PRESS     # Pressure in atm (default: 0.00101325)
--phi PHI               # Fugacity coefficient (default: 1.0)
```

#### Simulation Parameters
```bash
-n, --n-iterations N     # Number of hybrid iterations (default: 1000)
--write-interval N       # Snapshot save interval (default: 100)
```

#### Step Counts
```bash
--equil-steps N         # Empty framework equilibration (default: 100000)
--gcmc-steps N          # GCMC steps per iteration (default: 5000)
--nvt-steps N           # NVT relaxation steps (default: 10000)
--npt-steps N           # NPT MD steps (default: 50000)
```

## Output Files

### Directory Structure

```
T_87.3/
└── P_0.00101325/
    ├── current_config.data              # Current accepted configuration
    ├── acceptance_stats.txt             # Detailed iteration statistics
    ├── empty_framework_properties.txt   # U0 and V0 of empty framework
    ├── reference_config.data            # Latest NVT reference state
    │
    ├── config_iter_100.data             # Periodic snapshots
    ├── config_iter_200.data
    ├── config_iter_300.data
    │
    ├── log.gcmc.1.log                   # Individual step logs
    ├── log.gcmc.2.log
    ├── log.nvt.1.log
    ├── log.nvt.2.log
    ├── log.npt.1.log
    └── log.npt.2.log
```

## Restart Capability

pySIMAID automatically detects and resumes from `current_config.data`:

```bash
# Start simulation
./pysimaid.py -T 87.3 -P 0.00101325 -n 5000

# Interrupted at iteration 1234...
# [Ctrl+C or job timeout]

# Resume automatically (same command)
./pysimaid.py -T 87.3 -P 0.00101325 -n 5000
# Output: "Continuing from existing configuration"
# Will start from iteration 1235
```

**Note**: Statistics are appended to `acceptance_stats.txt`, so the full history is preserved.

## Thermodynamic Framework

### Free Energy in NPT Ensemble

The Gibbs free energy is computed as:

```
G = U + PV - TS
```

Where:
- **U** = potential energy
- **P** = pressure
- **V** = volume
- **T** = temperature
- **S** = configurational entropy

### Metropolis Acceptance Criterion
```
P_accept = min[1, exp(-beta DG)]

where beta = 1/(k_B T)
```

## Citation

If you use pySIMAID in your research, please cite:

```bibtex
@article{corrente2026pysimaid,
  title={A Thermodynamically Consistent Approach to Molecular Simulations of Adsorption-Induced Deformation and Structural Transitions in MOFs},
  author={Corrente, Nicholas J. and Chang, Kaelyn and Noor, Muhtasim and Neimark, Alexander V.},
  journal={ChemRxiv},
  year={2026},
  doi={10.26434/chemrxiv.15004165/v1},
  url={https://chemrxiv.org/doi/abs/10.26434/chemrxiv.15004165/v1}
}
```

## Todo list
- Make the atom types user input
- Mixtures
- Interface with open-source MC codes

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

