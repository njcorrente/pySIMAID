#!/usr/bin/env python3
"""
Hybrid MD/GCMC simulation controller with Metropolis acceptance criterion.

This script performs alternating GCMC and NPT MD simulations of gas adsorption in metal-organic frameworks. Each iteration consists of:
1. GCMC equilibration to sample particle insertions/deletions
2. NVT relaxation to remove hard overlaps on the full force field
3. NPT MD trajectory
4. Metropolis acceptance test based on free energy change

"""

import os
import sys
import subprocess
import shutil
import argparse
import random
import math
from pathlib import Path


class HybridSimulation:
    """Manages hybrid MD/GCMC simulation workflow."""
    
    def __init__(self, args):
        self.lammps_exec = args.lammps_exec
        self.nprocs = args.nprocs
        self.temperature = args.temperature
        self.pressure = args.pressure
        self.phi = args.phi
        self.n_iterations = args.n_iterations
        self.write_interval = args.write_interval
        
        # Simulation step counts
        self.equil_steps = args.equil_steps
        self.gcmc_steps = args.gcmc_steps
        self.nvt_steps = args.nvt_steps
        self.npt_steps = args.npt_steps
        
        # Physical constants
        self.kb = 0.001987204  # kcal/mol/K
        self.beta = 1.0 / (self.kb * self.temperature)
        
        # Directory setup
        self.home_dir = Path.cwd()
        self.run_dir = self.home_dir / f"T_{self.temperature}" / f"P_{self.pressure}"
        self.run_dir.mkdir(parents=True, exist_ok=True)
        
        # Statistics tracking
        self.n_accepted = 0
        self.n_rejected = 0
        self.stats_file = self.run_dir / "acceptance_stats.txt"
        
        # Check required input files
        self.check_input_files()
        
    def check_input_files(self):
        """Verify all required LAMMPS input files are present."""
        required_files = [
            "equilibrate_empty.in",
            "gcmc_step.in",
            "nvt_step.in",
            "npt_step.in",
            "paircoeffs.in",
            "data_zhang2013_SC"
        ]
        
        missing = [f for f in required_files if not Path(f).exists()]
        if missing:
            print("ERROR: Missing required input files:")
            for f in missing:
                print(f"  - {f}")
            sys.exit(1)
    
    def run_lammps(self, input_script, log_file, variables=None):
        """Execute LAMMPS with given input script and optional variables."""
        cmd = ["mpirun", "-np", str(self.nprocs), self.lammps_exec, 
               "-in", input_script, "-log", log_file]
        
        if variables:
            for key, value in variables.items():
                cmd.extend(["-var", key, str(value)])
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"ERROR: LAMMPS failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            sys.exit(1)
            
        return result
    
    def equilibrate_empty_framework(self):
        """Prepare initial structure and equilibrate empty framework at target T and P."""
        print(f"Equilibrating empty framework at T={self.temperature} K, "
              f"P={self.pressure} atm for {self.equil_steps} steps...")
        
        variables = {
            'Tsim': self.temperature,
            'p': self.pressure,
            'equil_steps': self.equil_steps
        }
        
        self.run_lammps("equilibrate_empty.in", "equilibrate_empty.log", variables)
        
        # Read and store empty framework properties
        with open("empty_framework_properties.txt") as f:
            v0, u0 = f.read().strip().split()
        print(f"Empty framework equilibrated: V0 = {v0} A^3, U0 = {u0} kcal/mol")
        
        # Copy equilibrated structure to run directory
        shutil.copy("empty_equilibrated.data", self.run_dir / "current_config.data")
        shutil.copy("empty_framework_properties.txt", self.run_dir)
        
    def check_restart(self):
        """Check if we're continuing from a previous run."""
        current_config = self.run_dir / "current_config.data"
        
        if current_config.exists():
            print("Continuing from existing configuration")
            return True
        else:
            print("Starting fresh simulation")
            return False
            
    def extract_final_state(self, state_file):
        """Parse final thermodynamic state from LAMMPS output."""
        with open(state_file) as f:
            line = f.readlines()[-1].strip()
            values = line.split()
            return {
                'u': float(values[0]),
                'v': float(values[1]),
                'p': float(values[2]),
                'g': float(values[3])
            }
    
    def metropolis_accept(self, delta_g):
        """Apply Metropolis criterion for move acceptance."""
        rand = random.random()
        
        if delta_g <= 0:
            acc_prob = 1.0
        else:
            acc_prob = math.exp(-self.beta * delta_g)
        
        accept = rand < acc_prob
        
        return accept, acc_prob, rand
    
    def run_gcmc_step(self, iteration):
        """Execute GCMC equilibration step."""
        print(f"[GCMC] Starting GCMC equilibration...")
        
        # Copy current accepted config as GCMC input
        shutil.copy(self.run_dir / "current_config.data", 
                   self.run_dir / "gcmc_initial.data")
        
        variables = {
            'Tsim': self.temperature,
            'p': self.pressure,
            'phi': self.phi,
            'rundir': str(self.run_dir),
            'iter': iteration,
            'gcmc_steps': self.gcmc_steps
        }
        
        self.run_lammps("gcmc_step.in", 
                       str(self.run_dir / f"log.gcmc.{iteration}.log"),
                       variables)
        
        state = self.extract_final_state(self.run_dir / "gcmc_final_state.txt")
        print(f"[GCMC] Final: U={state['u']:.2f} V={state['v']:.2f} G={state['g']:.2f}")
        
        return state
    
    def run_nvt_step(self, iteration):
        """Execute NVT relaxation to remove hard overlaps."""
        print(f"[NVT] Relaxing GCMC output on full force field...")
        
        variables = {
            'Tsim': self.temperature,
            'rundir': str(self.run_dir),
            'iter': iteration,
            'nvt_steps': self.nvt_steps
        }
        
        self.run_lammps("nvt_step.in",
                       str(self.run_dir / f"log.nvt.{iteration}.log"),
                       variables)
        
        # NVT output becomes both NPT input and reference state
        shutil.copy(self.run_dir / "nvt_final.data",
                   self.run_dir / "npt_initial.data")
        shutil.copy(self.run_dir / "nvt_final.data",
                   self.run_dir / "reference_config.data")
        
        # Read reference state properties
        with open(self.run_dir / "reference_u.txt") as f:
            u_ref = float(f.read().strip())
        with open(self.run_dir / "reference_v.txt") as f:
            v_ref = float(f.read().strip())
        with open(self.run_dir / "reference_n.txt") as f:
            n_ref = int(f.read().strip())
        
        print(f"[REFERENCE] U_ref = {u_ref:.2f} kcal/mol, "
              f"V_ref = {v_ref:.2f} A^3, N_ref = {n_ref}")
        
        return {'u': u_ref, 'v': v_ref, 'n': n_ref}
    
    def run_npt_step(self, iteration):
        """Execute NPT MD simulation."""
        print(f"[NPT] Starting NPT MD simulation...")
        
        variables = {
            'Tsim': self.temperature,
            'p': self.pressure,
            'rundir': str(self.run_dir),
            'iter': iteration,
            'kb': self.kb,
            'npt_steps': self.npt_steps
        }
        
        self.run_lammps("npt_step.in",
                       str(self.run_dir / f"log.npt.{iteration}.log"),
                       variables)
        
        state = self.extract_final_state(self.run_dir / "npt_final_state.txt")
        print(f"[NPT] Final: U={state['u']:.2f} V={state['v']:.2f} G={state['g']:.2f}")
        
        return state
    
    def apply_acceptance(self, iteration, g_gcmc, g_npt):
        """Apply Metropolis acceptance criterion and update configuration."""
        print(f"[METROPOLIS] Applying acceptance criterion...")
        print(f"[METROPOLIS] Comparing NPT final vs NVT output (reference state)")
        
        accept, acc_prob, rand_num = self.metropolis_accept(g_npt)
        
        print(f"[METROPOLIS] ΔG_GCMC = {g_gcmc:.4f} kcal/mol (informational)")
        print(f"[METROPOLIS] ΔG_NPT = {g_npt:.4f} kcal/mol (used for acceptance)")
        print(f"[METROPOLIS] Acceptance probability = {acc_prob:.4f}")
        print(f"[METROPOLIS] Random number = {rand_num:.4f}")
        
        if accept:
            print("[METROPOLIS] ✓ ACCEPTED - Updating to NPT final configuration")
            shutil.copy(self.run_dir / "npt_final.data",
                       self.run_dir / "current_config.data")
            self.n_accepted += 1
            status = "ACCEPTED"
        else:
            print("[METROPOLIS] ✗ REJECTED - Keeping NVT-relaxed output (state B)")
            # GCMC + NVT moves are kept; only NPT trajectory is rejected
            shutil.copy(self.run_dir / "nvt_final.data",
                       self.run_dir / "current_config.data")
            self.n_rejected += 1
            status = "REJECTED"
        
        # Log statistics
        with open(self.stats_file, "a") as f:
            f.write(f"{iteration} {int(accept)} {g_gcmc:.6f} {g_npt:.6f} "
                   f"{acc_prob:.6f} {rand_num:.6f}\n")
        
        return status
    
    def write_periodic_output(self, iteration):
        """Save configuration snapshot."""
        print(f"[OUTPUT] Writing configuration snapshot...")
        shutil.copy(self.run_dir / "current_config.data",
                   self.run_dir / f"config_iter_{iteration}.data")
    
    def print_statistics(self, status):
        """Print running acceptance statistics."""
        total = self.n_accepted + self.n_rejected
        acc_rate = 100.0 * self.n_accepted / total if total > 0 else 0
        
        print("\nRunning Statistics:")
        print(f"  Status: {status}")
        print(f"  Accepted: {self.n_accepted}")
        print(f"  Rejected: {self.n_rejected}")
        print(f"  Acceptance Rate: {acc_rate:.2f}%")
        print("=" * 60)
    
    def run(self):
        """Execute the full simulation workflow."""
        print("=" * 60)
        print("Hybrid MD/GCMC Simulation")
        print(f"T={self.temperature} K, P={self.pressure} atm, phi={self.phi}")
        print(f"N_ITERATIONS={self.n_iterations}")
        print(f"Run lengths: EQUIL={self.equil_steps} GCMC={self.gcmc_steps} "
              f"NVT={self.nvt_steps} NPT={self.npt_steps}")
        print("=" * 60)
        
        # Check if continuing from previous run
        is_restart = self.check_restart()
        
        if not is_restart:
            # Fresh start: prepare and equilibrate empty framework
            self.equilibrate_empty_framework()
            
            # Initialize statistics file
            with open(self.stats_file, "w") as f:
                f.write("# Iter Accept G_GCMC G_NPT AccProb RandNum\n")
        
        # Main iteration loop
        for iteration in range(1, self.n_iterations + 1):
            print(f"\n{'=' * 20} Iteration {iteration}/{self.n_iterations} {'=' * 20}")
            
            # Step 1: GCMC equilibration
            gcmc_state = self.run_gcmc_step(iteration)
            
            # Step 1.5: NVT relaxation (produces reference state for NPT)
            nvt_state = self.run_nvt_step(iteration)
            
            # Step 2: NPT MD trajectory
            npt_state = self.run_npt_step(iteration)
            
            # Step 3: Metropolis acceptance test
            status = self.apply_acceptance(iteration, gcmc_state['g'], npt_state['g'])
            
            # Step 4: Periodic output
            if iteration % self.write_interval == 0:
                self.write_periodic_output(iteration)
            
            # Print statistics
            self.print_statistics(status)
        
        # Final summary
        total = self.n_accepted + self.n_rejected
        final_acc_rate = 100.0 * self.n_accepted / total if total > 0 else 0
        
        print("\n" + "=" * 60)
        print("Simulation Complete!")
        print("=" * 60)
        print(f"Total iterations: {self.n_iterations}")
        print(f"Accepted: {self.n_accepted}")
        print(f"Rejected: {self.n_rejected}")
        print(f"Final acceptance rate: {final_acc_rate:.2f}%")
        print("=" * 60)
        print(f"\nResults saved to: {self.run_dir}")
        print("Log files: log.gcmc.*.log, log.nvt.*.log, log.npt.*.log")
        print("Acceptance stats: acceptance_stats.txt")


def main():
    parser = argparse.ArgumentParser(
        description="Hybrid MD/GCMC simulation with Metropolis acceptance",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # LAMMPS execution parameters
    parser.add_argument("--lammps-exec", default="lmp_mpi",
                       help="Path to LAMMPS executable")
    parser.add_argument("--nprocs", type=int, default=1,
                       help="Number of MPI processes")
    
    # Thermodynamic conditions
    parser.add_argument("-T", "--temperature", type=float, default=87.3,
                       help="Temperature (K)")
    parser.add_argument("-P", "--pressure", type=float, default=0.00101325,
                       help="Pressure (atm)")
    parser.add_argument("--phi", type=float, default=1.0,
                       help="Fugacity coefficient")
    
    # Simulation parameters
    parser.add_argument("-n", "--n-iterations", type=int, default=1000,
                       help="Number of hybrid iterations")
    parser.add_argument("--write-interval", type=int, default=100,
                       help="Snapshot save interval")
    
    # Step counts
    parser.add_argument("--equil-steps", type=int, default=100000,
                       help="Empty framework equilibration steps")
    parser.add_argument("--gcmc-steps", type=int, default=5000,
                       help="GCMC steps per iteration")
    parser.add_argument("--nvt-steps", type=int, default=10000,
                       help="NVT relaxation steps per iteration")
    parser.add_argument("--npt-steps", type=int, default=50000,
                       help="NPT MD steps per iteration")
    
    args = parser.parse_args()
    
    # Run simulation
    sim = HybridSimulation(args)
    sim.run()


if __name__ == "__main__":
    main()