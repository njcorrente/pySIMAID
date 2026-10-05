# pySIMAID — Simulator for Adsorption-Induced Deformation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![LAMMPS](https://img.shields.io/badge/LAMMPS-required-red.svg)](https://www.lammps.org/)

A Python-controlled workflow for simulating gas adsorption and structural deformation in metal-organic frameworks (MOFs) using a thermodynamically rigorous hybrid MD/GCMC approach.

## Overview

**pySIMAID** (Python Simulator for Adsorption-Induced Deformation) captures the coupled phenomena of gas adsorption and framework flexibility in nanoporous materials. The method combines:

- **GCMC** for efficient sampling of particle insertions and deletions at fixed chemical potential
- **NVT relaxation** for adsorbate relaxation and velocity assignment
- **NPT molecular dynamics** for structural relaxation and volume fluctuations
- **Metropolis acceptance** based on free-energy comparisons
- **Post-processing tools** for extracting adsorption and deformation results from LAMMPS output

This approach enables the prediction of adsorption isotherms and framework deformation in flexible MOFs.

## Key Features

- **Thermodynamically consistent** acceptance criterion
- **Automatic restart** from existing checkpoints
- **Flexible configuration** through command-line arguments
- **Detailed statistics** tracking and logging
- **Modular LAMMPS inputs** for force-field customization
- **MPI-parallel** execution
- **Optional Metropolis filtering**, allowing every NPT trajectory to be accepted
- **Automated result extraction** from hybrid simulation directories
- **Per-cycle and bootstrap output** for statistical analysis

## Scientific Motivation

Traditional GCMC simulations generally assume rigid frameworks, while purely MD-based approaches may require unfeasibly large systems and simulation times. pySIMAID bridges this gap by:

1. Using GCMC to efficiently sample adsorption and desorption
2. Allowing the framework to relax and deform through molecular dynamics
3. Maintaining detailed balance through Metropolis acceptance
4. Extracting adsorption and structural properties across hybrid cycles

This workflow can capture adsorption-induced deformation phenomena such as:

- Framework breathing and swelling
- Gate-opening transitions
- Cooperative adsorption effects
- Pressure-dependent structural changes

## Requirements

### Software Dependencies

- **Python 3.7+**
- **LAMMPS** (2020 or later) compiled with:
  - `MC` package for GCMC
  - `KSPACE` package for long-range electrostatics
  - `MOLECULE` package for molecular systems
  - MPI support
- An **MPI implementation**, such as:
  - OpenMPI
  - MPICH
  - Intel MPI

Install any Python dependencies listed in `requirements.txt`:

```bash
python -m pip install -r requirements.txt
```

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/pySIMAID.git
cd pySIMAID

# Install Python dependencies
python -m pip install -r requirements.txt

# Verify the LAMMPS installation
mpirun -np 1 lmp_mpi -help

# Make the main script executable
chmod +x pysimaid.py

# Test the simulation driver
./pysimaid.py --help

# Test the result extractor
python tools/extract_hybrid_results.py --help
```

## Quick Start

Run a simulation at 87.3 K and 0.00101325 atm:

```bash
./pysimaid.py -T 87.3 -P 0.00101325 -n 1000
```

Run without the Metropolis acceptance criterion:

```bash
./pysimaid.py \
    -T 87.3 \
    -P 0.00101325 \
    -n 1000 \
    --no-metropolis
```

Monitor the simulation:

```bash
tail -f T_87.3/P_0.00101325/acceptance_stats.txt
```

Count accepted hybrid cycles:

```bash
grep -c "1$" T_87.3/P_0.00101325/acceptance_stats.txt
```

Extract adsorption and deformation results:

```bash
python tools/extract_hybrid_results.py \
    T_87.3 \
    structures/your_framework.data \
    --temperature 87.3 \
    --accepted-only \
    --per-cycle-csv per_cycle.csv \
    --bootstrap-csv bootstrap.csv \
    --output isotherm.csv
```

## Repository Structure

```text
pySIMAID/
├── pysimaid.py                       # Main simulation control script
│
├── tools/
│   ├── extract_hybrid_results.py     # Result extraction utility
│   └── LammpsLogParser.py            # LAMMPS log parser
│
├── lammps_inputs/
│   ├── equilibrate_empty.in          # Prepare and equilibrate the empty framework
│   ├── gcmc_step.in                  # GCMC equilibration
│   ├── nvt_step.in                   # NVT relaxation and reference-state generation
│   ├── npt_step.in                   # NPT molecular-dynamics trajectory
│   └── paircoeffs.in                 # Force-field parameters
│
├── structures/
│   └── *.data                        # Initial MOF structures
│
├── README.md
├── LICENSE
└── requirements.txt
```

The extractor and its parser are kept together in `tools/` so that the following import works without modifying `PYTHONPATH`:

```python
from LammpsLogParser import LammpsLogParser
```

Run the extractor from the repository root with:

```bash
python tools/extract_hybrid_results.py ...
```

## Simulation Usage

### Basic Command

```bash
./pysimaid.py [OPTIONS]
```

### LAMMPS Execution

```text
--lammps-exec PATH       Path to the LAMMPS executable
                         Default: lmp_mpi

--nprocs N               Number of MPI processes
                         Default: 16
```

### Thermodynamic Conditions

```text
-T, --temperature TEMP   Temperature in kelvin
                         Default: 87.3

-P, --pressure PRESS     Pressure in atm
                         Default: 0.00101325

--phi PHI                Fugacity coefficient
                         Default: 1.0
```

### Simulation Parameters

```text
-n, --n-iterations N     Number of hybrid iterations
                         Default: 1000

--write-interval N       Snapshot save interval
                         Default: 100

--no-metropolis          Disable Metropolis filtering and accept every
                         NPT result
```

### Step Counts

```text
--equil-steps N          Empty-framework equilibration steps
                         Default: 100000

--gcmc-steps N           GCMC steps per hybrid iteration
                         Default: 5000

--nvt-steps N            NVT relaxation steps per hybrid iteration
                         Default: 10000

--npt-steps N            NPT MD steps per hybrid iteration
                         Default: 50000
```

For the complete and authoritative option list, run:

```bash
./pysimaid.py --help
```

## Simulation Output

### Directory Structure

A typical simulation produces the following directory hierarchy:

```text
T_87.3/
└── P_0.00101325/
    ├── current_config.data              # Current accepted configuration
    ├── acceptance_stats.txt             # Hybrid-cycle statistics
    ├── empty_framework_properties.txt   # Empty-framework U0 and V0
    ├── reference_config.data            # Most recent NVT reference state
    │
    ├── config_iter_100.data             # Periodic snapshots
    ├── config_iter_200.data
    ├── config_iter_300.data
    │
    ├── log.gcmc.1.log                   # GCMC logs
    ├── log.gcmc.2.log
    ├── log.nvt.1.log                    # NVT logs
    ├── log.nvt.2.log
    ├── log.npt.1.log                    # NPT logs
    └── log.npt.2.log
```

## Restart Capability

pySIMAID automatically detects and resumes from `current_config.data`.

Start a simulation:

```bash
./pysimaid.py -T 87.3 -P 0.00101325 -n 5000
```

If the simulation is interrupted, rerun the same command:

```bash
./pysimaid.py -T 87.3 -P 0.00101325 -n 5000
```

Example restart message:

```text
Continuing from existing configuration
```

If the previous run ended at iteration 1234, the restarted simulation continues from iteration 1235.

Statistics are appended to `acceptance_stats.txt`, preserving the simulation history across restarts.

## Result Extraction

The post-processing utility is located at:

```text
tools/extract_hybrid_results.py
```

Its associated LAMMPS parser is located at:

```text
tools/LammpsLogParser.py
```

### Basic Extraction

```bash
python tools/extract_hybrid_results.py \
    PATH_TO_RESULTS \
    structures/your_framework.data \
    --temperature 87.3 \
    --output isotherm.csv
```

### Accepted Configurations Only

Use `--accepted-only` to restrict the analysis to accepted hybrid configurations:

```bash
python tools/extract_hybrid_results.py \
    PATH_TO_RESULTS \
    structures/your_framework.data \
    --temperature 87.3 \
    --accepted-only \
    --output isotherm.csv
```

### Per-Cycle and Bootstrap Results

```bash
python tools/extract_hybrid_results.py \
    PATH_TO_RESULTS \
    structures/your_framework.data \
    --temperature 87.3 \
    --accepted-only \
    --per-cycle-csv per_cycle.csv \
    --bootstrap-csv bootstrap.csv \
    --output isotherm.csv
```

The output files serve different purposes:

- `isotherm.csv` contains aggregated adsorption and deformation results.
- `per_cycle.csv` contains values extracted for individual hybrid cycles.
- `bootstrap.csv` contains bootstrap-based statistical estimates.

For the complete list of extractor arguments, run:

```bash
python tools/extract_hybrid_results.py --help
```

## Thermodynamic Framework

### Free Energy in the NPT Ensemble

The Gibbs free energy is represented as:

```text
G = U + PV - TS
```

where:

- `U` is the potential energy
- `P` is the pressure
- `V` is the volume
- `T` is the temperature
- `S` is the configurational entropy

### Metropolis Acceptance Criterion

By default, the final configuration generated by each NPT trajectory is accepted according to:

```text
P_accept = min[1, exp(-βΔG)]
```

where:

```text
β = 1 / (k_B T)
```

If the NPT move is accepted, `npt_final.data` becomes the configuration used by the next hybrid iteration. If it is rejected, the NVT-relaxed reference configuration is retained.

### Disabling Metropolis Filtering

The Metropolis criterion can be disabled with:

```bash
./pysimaid.py --no-metropolis
```

In this mode:

```text
P_accept = 1
```

Every NPT final configuration is therefore accepted. The GCMC, NVT, and NPT stages are otherwise unchanged, and acceptance statistics continue to be written to `acceptance_stats.txt`.

The default behavior is unchanged when `--no-metropolis` is not supplied.

## Citation

If you use pySIMAID in your research, please cite:

```bibtex
@article{corrente2026pysimaid,
  title   = {A Thermodynamically Consistent Approach to Molecular Simulations of Adsorption-Induced Deformation and Structural Transitions in MOFs},
  author  = {Corrente, Nicholas J. and Chang, Kaelyn and Noor, Muhtasim and Neimark, Alexander V.},
  journal = {ChemRxiv},
  year    = {2026},
  doi     = {10.26434/chemrxiv.15004165/v1},
  url     = {https://chemrxiv.org/doi/abs/10.26434/chemrxiv.15004165/v1}
}
```

## Todo

- Add mixture examples
- Add an interface to open-source Monte Carlo codes
- Add example post-processing datasets
- Add automated tests for the simulation and extraction workflows

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
