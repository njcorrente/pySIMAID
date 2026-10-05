#!/usr/bin/env python3
"""Parse thermo sections from LAMMPS log files."""

from pathlib import Path

import numpy as np


class LammpsLogParser:
    """Parse LAMMPS thermo output and calculate adsorbate loading."""

    INT_COLUMNS = {"Step", "Atoms"}

    def __init__(self, log_file_path=None):
        self.log_file_path = None
        self.results = {}
        self.structure_atoms = None
        self.framework_mass_g = None

        if log_file_path is not None:
            self.parse_file(log_file_path)

    def parse_file(self, log_file_path=None):
        if log_file_path is not None:
            self.log_file_path = Path(log_file_path)

        if self.log_file_path is None:
            raise ValueError("No LAMMPS log file was provided")

        self.results = {}
        columns = None
        rows = []
        section_number = 0
        in_thermo = False

        with self.log_file_path.open("r", errors="replace") as handle:
            for line in handle:
                stripped = line.strip()

                if self._is_thermo_header(stripped):
                    if columns is not None and rows:
                        section_number += 1
                        self.results[f"thermo_{section_number}"] = rows

                    columns = stripped.split()
                    rows = []
                    in_thermo = True
                    continue

                if in_thermo and stripped.startswith("Loop time of"):
                    in_thermo = False
                    continue

                if in_thermo and columns:
                    row = self._parse_thermo_row(stripped, columns)
                    if row is not None:
                        rows.append(row)

        if columns is not None and rows:
            section_number += 1
            self.results[f"thermo_{section_number}"] = rows

        return self.results

    def set_framework_datafile(self, data_file_path):
        """Read framework atom counts and masses from a LAMMPS data file."""
        masses, atom_counts = self._parse_data_file(data_file_path)

        missing_masses = sorted(set(atom_counts) - set(masses))
        if missing_masses:
            missing = ", ".join(str(value) for value in missing_masses)
            raise ValueError(
                f"Missing masses for atom type(s) {missing} in "
                f"{data_file_path}"
            )

        mass_amu = sum(
            masses[atom_type] * count
            for atom_type, count in atom_counts.items()
        )

        if mass_amu <= 0:
            raise ValueError(
                f"Could not determine framework mass from {data_file_path}"
            )

        self.structure_atoms = sum(atom_counts.values())
        self.framework_mass_g = mass_amu * 1.66053906660e-24
        self._add_adsorption_columns()

    def get_chunk_names(self):
        return list(self.results)

    def get_chunk_as_arrays(self, chunk_name):
        rows = self.results.get(chunk_name, [])
        if not rows:
            return {}

        keys = rows[0].keys()
        return {
            key: np.asarray([row[key] for row in rows])
            for key in keys
        }

    @staticmethod
    def _is_thermo_header(line):
        parts = line.split()
        if len(parts) < 2 or parts[0] != "Step":
            return False

        try:
            float(parts[1])
        except ValueError:
            return True

        return False

    def _parse_thermo_row(self, line, columns):
        parts = line.split()
        if len(parts) != len(columns):
            return None

        row = {}
        try:
            for column, value in zip(columns, parts):
                if column in self.INT_COLUMNS:
                    row[column] = int(value)
                else:
                    row[column] = float(value)
        except ValueError:
            return None

        return row

    def _add_adsorption_columns(self):
        if self.structure_atoms is None or self.framework_mass_g is None:
            return

        avogadro = 6.02214076e23

        for rows in self.results.values():
            for row in rows:
                if "Atoms" not in row:
                    continue

                count = row["Atoms"] - self.structure_atoms
                row["AdsorptionCount"] = count
                row["AdsorptionMmolPerGram"] = (
                    count * 1000.0
                    / (avogadro * self.framework_mass_g)
                )

    @staticmethod
    def _parse_data_file(data_file_path):
        """Read Masses and Atoms sections from an atom_style full data file."""
        masses = {}
        atom_counts = {}
        section = None

        with Path(data_file_path).open("r") as handle:
            for line in handle:
                content = line.split("#", 1)[0].strip()

                if not content:
                    continue

                if content == "Masses":
                    section = "Masses"
                    continue

                if content == "Atoms":
                    section = "Atoms"
                    continue

                parts = content.split()

                try:
                    int(parts[0])
                except (ValueError, IndexError):
                    section = None
                    continue

                if section == "Masses" and len(parts) >= 2:
                    masses[int(parts[0])] = float(parts[1])

                elif section == "Atoms" and len(parts) >= 3:
                    # atom_style full: atom-ID molecule-ID atom-type ...
                    atom_type = int(parts[2])
                    atom_counts[atom_type] = (
                        atom_counts.get(atom_type, 0) + 1
                    )

        if not atom_counts:
            raise ValueError(
                f"No atoms were found in the Atoms section of {data_file_path}"
            )

        return masses, atom_counts
