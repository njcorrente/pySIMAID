#!/usr/bin/env python3
"""Extract adsorption and volume statistics from pySIMAID simulations."""

import argparse
import csv
import math
import re
from pathlib import Path

import numpy as np

from LammpsLogParser import LammpsLogParser


KB_SI = 1.380649e-23
ANGSTROM3_TO_M3 = 1.0e-30
ITERATION_PATTERN = re.compile(r"\.(\d+)\.log$")


def iteration_number(path):
    match = ITERATION_PATTERN.search(Path(path).name)
    if match is None:
        raise ValueError(f"Cannot determine iteration number from {path}")
    return int(match.group(1))


def find_pressure_directories(root_dir, temperature=None):
    root = Path(root_dir)
    directories = []

    for pressure_dir in root.rglob("P_*"):
        if not pressure_dir.is_dir():
            continue

        temperature_dir = pressure_dir.parent
        if not temperature_dir.name.startswith("T_"):
            continue

        try:
            found_temperature = float(temperature_dir.name[2:])
            pressure = float(pressure_dir.name[2:])
        except ValueError:
            continue

        if temperature is not None and not math.isclose(
            found_temperature,
            float(temperature),
            rel_tol=1.0e-9,
            abs_tol=1.0e-12,
        ):
            continue

        directories.append((pressure, pressure_dir))

    directories.sort(key=lambda item: item[0])
    return directories


def find_iteration_pairs(pressure_dir):
    gcmc_logs = {
        iteration_number(path): path
        for path in pressure_dir.glob("log.gcmc.*.log")
    }
    npt_logs = {
        iteration_number(path): path
        for path in pressure_dir.glob("log.npt.*.log")
    }

    common = sorted(set(gcmc_logs) & set(npt_logs))
    missing_gcmc = sorted(set(npt_logs) - set(gcmc_logs))
    missing_npt = sorted(set(gcmc_logs) - set(npt_logs))

    if missing_gcmc:
        print(
            f"  Warning: missing GCMC logs for iterations "
            f"{format_iteration_list(missing_gcmc)}"
        )

    if missing_npt:
        print(
            f"  Warning: missing NPT logs for iterations "
            f"{format_iteration_list(missing_npt)}"
        )

    return [
        (iteration, gcmc_logs[iteration], npt_logs[iteration])
        for iteration in common
    ]


def format_iteration_list(iterations, limit=10):
    if len(iterations) <= limit:
        return ", ".join(str(value) for value in iterations)

    head = ", ".join(str(value) for value in iterations[:limit])
    return f"{head}, ... ({len(iterations)} total)"


def load_acceptance_stats(pressure_dir):
    path = pressure_dir / "acceptance_stats.txt"
    if not path.exists():
        return {}

    records = {}

    with path.open("r") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            parts = stripped.split()
            if len(parts) < 2:
                continue

            try:
                iteration = int(parts[0])
                accepted = bool(int(parts[1]))
            except ValueError:
                print(
                    f"  Warning: invalid acceptance record at "
                    f"{path}:{line_number}"
                )
                continue

            if iteration in records:
                print(
                    f"  Warning: duplicate acceptance record for iteration "
                    f"{iteration}; using the last value"
                )

            records[iteration] = accepted

    return records


def load_production_section(log_path, framework_datafile):
    parser = LammpsLogParser(log_path)
    parser.set_framework_datafile(framework_datafile)

    chunks = parser.get_chunk_names()
    if not chunks:
        raise ValueError(f"No thermo data found in {log_path}")

    # pySIMAID writes a run-0 section followed by the production run.
    return parser.get_chunk_as_arrays(chunks[-1])


def find_column(arrays, names):
    for name in names:
        if name in arrays:
            return arrays[name]
    return None


def select_samples(values, equilibration_fraction, downsample):
    values = np.asarray(values, dtype=float)

    if values.size == 0:
        return values

    start = int(values.size * equilibration_fraction)
    return values[start::downsample]


def volume_moments(volume_angstrom3):
    values = np.asarray(volume_angstrom3, dtype=float)

    if values.size == 0:
        return np.nan, np.nan, np.nan

    mean = float(np.mean(values))
    mean_square = float(np.mean(values * values))
    variance = max(mean_square - mean * mean, 0.0)

    return mean, mean_square, variance


def compressibility(volume_angstrom3, temperature):
    mean, _, variance = volume_moments(volume_angstrom3)

    if not np.isfinite(mean) or mean <= 0 or temperature <= 0:
        return np.nan

    mean_m3 = mean * ANGSTROM3_TO_M3
    variance_m6 = variance * ANGSTROM3_TO_M3**2

    return variance_m6 / (KB_SI * temperature * mean_m3)


def integrated_autocorrelation_time(values):
    values = np.asarray(values, dtype=float)
    size = values.size

    if size < 10:
        return 1.0

    centered = values - np.mean(values)
    variance = np.dot(centered, centered)

    if variance <= 0:
        return 1.0

    fft_size = 1 << (2 * size - 1).bit_length()
    transformed = np.fft.rfft(centered, fft_size)
    correlation = np.fft.irfft(
        transformed * np.conjugate(transformed),
        fft_size,
    )[:size]

    correlation /= correlation[0]

    negative = np.flatnonzero(correlation[1:] <= 0)
    cutoff = negative[0] + 1 if negative.size else size

    tau = 1.0 + 2.0 * np.sum(correlation[1:cutoff])
    return float(max(1.0, min(tau, size / 2.0)))


def bootstrap_compressibility(
    volume_chunks_angstrom3,
    temperature,
    n_bootstrap,
    rng,
):
    chunks = [
        np.asarray(chunk, dtype=float) * ANGSTROM3_TO_M3
        for chunk in volume_chunks_angstrom3
        if chunk is not None and len(chunk) > 0
    ]

    n_chunks = len(chunks)
    if n_chunks == 0:
        return {
            "kappa": np.nan,
            "kappa_std": np.nan,
            "tau_iterations": np.nan,
            "block_size": 0,
            "effective_chunks": np.nan,
        }

    all_values = np.concatenate(chunks)
    mean = np.mean(all_values)

    if mean <= 0:
        kappa = np.nan
    else:
        kappa = np.var(all_values) / (KB_SI * temperature * mean)

    if n_chunks < 2 or n_bootstrap < 2:
        return {
            "kappa": float(kappa),
            "kappa_std": np.nan,
            "tau_iterations": 1.0,
            "block_size": 1,
            "effective_chunks": float(n_chunks),
        }

    chunk_means = np.asarray([np.mean(chunk) for chunk in chunks])
    tau = integrated_autocorrelation_time(chunk_means)
    block_size = min(n_chunks, max(1, math.ceil(2.0 * tau)))
    effective_chunks = n_chunks / (2.0 * tau)

    n_starts = n_chunks - block_size + 1
    n_blocks = math.ceil(n_chunks / block_size)
    estimates = []

    for _ in range(n_bootstrap):
        sampled_chunks = []

        for start in rng.integers(0, n_starts, size=n_blocks):
            sampled_chunks.extend(chunks[start:start + block_size])

        sampled = np.concatenate(sampled_chunks[:n_chunks])
        sampled_mean = np.mean(sampled)

        if sampled_mean > 0:
            estimates.append(
                np.var(sampled)
                / (KB_SI * temperature * sampled_mean)
            )

    if len(estimates) > 1:
        kappa_std = float(np.std(estimates, ddof=1))
    else:
        kappa_std = np.nan

    return {
        "kappa": float(kappa),
        "kappa_std": kappa_std,
        "tau_iterations": float(tau),
        "block_size": int(block_size),
        "effective_chunks": float(effective_chunks),
    }


def tail_fraction_for_pressure(
    pressure,
    tail_fraction,
    tail_fraction_low,
    pressure_switch,
):
    if (
        tail_fraction_low is not None
        and pressure_switch is not None
        and pressure < pressure_switch
    ):
        return tail_fraction_low

    return tail_fraction


def process_pressure(
    pressure,
    pressure_dir,
    framework_datafile,
    temperature,
    equilibration_fraction,
    tail_fraction,
    tail_fraction_low,
    pressure_switch,
    downsample,
    accepted_only,
    n_bootstrap,
    random_seed,
):
    print(f"Processing {pressure_dir}")

    pairs = find_iteration_pairs(pressure_dir)
    if not pairs:
        print("  Warning: no paired GCMC/NPT logs found")
        return None

    acceptance = load_acceptance_stats(pressure_dir)

    if accepted_only:
        if not acceptance:
            print(
                "  Warning: --accepted-only was requested, but "
                "acceptance_stats.txt was not found or was empty"
            )
            return None

        missing = [
            iteration
            for iteration, _, _ in pairs
            if iteration not in acceptance
        ]
        if missing:
            print(
                f"  Warning: no acceptance result for iterations "
                f"{format_iteration_list(missing)}; excluding them"
            )

        pairs = [
            pair
            for pair in pairs
            if acceptance.get(pair[0], False)
        ]

        print(f"  Accepted iteration pairs: {len(pairs)}")

    if not pairs:
        print("  Warning: no iterations remain after acceptance filtering")
        return None

    fraction = tail_fraction_for_pressure(
        pressure,
        tail_fraction,
        tail_fraction_low,
        pressure_switch,
    )
    n_use = max(1, math.ceil(len(pairs) * fraction))
    pairs = pairs[-n_use:]

    print(f"  Using the final {n_use} paired iterations")

    adsorption_chunks = []
    volume_chunks = []
    per_cycle = []

    for iteration, gcmc_log, npt_log in pairs:
        try:
            gcmc = load_production_section(
                gcmc_log,
                framework_datafile,
            )
            npt = load_production_section(
                npt_log,
                framework_datafile,
            )

            adsorption_full = find_column(
                gcmc,
                ("AdsorptionMmolPerGram",),
            )
            volume_full = find_column(
                npt,
                ("Volume", "v_volume"),
            )

            if adsorption_full is None:
                raise ValueError(
                    "AdsorptionMmolPerGram is absent from the GCMC log"
                )

            if volume_full is None:
                raise ValueError(
                    "Volume or v_volume is absent from the NPT log"
                )

            adsorption = select_samples(
                adsorption_full,
                equilibration_fraction,
                downsample,
            )
            volume = select_samples(
                volume_full,
                equilibration_fraction,
                downsample,
            )

            if adsorption.size == 0 or volume.size == 0:
                raise ValueError(
                    "No samples remain after equilibration trimming"
                )

            # The NPT stage starts from the final GCMC loading.
            n_alpha = float(adsorption_full[-1])

            v_mean, v_mean_square, v_variance = volume_moments(volume)
            kappa_tn = compressibility(volume, temperature)

            adsorption_chunks.append(adsorption)
            volume_chunks.append(volume)

            per_cycle.append({
                "pressure_atm": pressure,
                "iteration": iteration,
                "accepted": acceptance.get(iteration),
                "N_alpha_mmol_g": n_alpha,
                "volume_mean_A3": v_mean,
                "volume_mean_square_A6": v_mean_square,
                "volume_variance_A6": v_variance,
                "kappa_TN_Pa_inverse": kappa_tn,
                "n_volume_samples": int(volume.size),
            })

        except (OSError, ValueError) as error:
            print(f"  Warning: iteration {iteration}: {error}")

    if not adsorption_chunks or not volume_chunks:
        print("  Warning: no usable data were extracted")
        return None

    adsorption_all = np.concatenate(adsorption_chunks)
    volume_all = np.concatenate(volume_chunks)

    bootstrap = bootstrap_compressibility(
        volume_chunks,
        temperature,
        n_bootstrap,
        np.random.default_rng(random_seed),
    )

    result = {
        "pressure_atm": pressure,
        "adsorption_mean_mmol_g": float(np.mean(adsorption_all)),
        "adsorption_std_mmol_g": (
            float(np.std(adsorption_all, ddof=1))
            if adsorption_all.size > 1
            else 0.0
        ),
        "volume_mean_A3": float(np.mean(volume_all)),
        "volume_std_A3": (
            float(np.std(volume_all, ddof=1))
            if volume_all.size > 1
            else 0.0
        ),
        "volume_variance_A6": float(np.var(volume_all)),
        "kappa_Tmu_Pa_inverse": bootstrap["kappa"],
        "kappa_Tmu_std_Pa_inverse": bootstrap["kappa_std"],
        "n_iterations": len(volume_chunks),
        "n_adsorption_samples": int(adsorption_all.size),
        "n_volume_samples": int(volume_all.size),
    }

    diagnostics = {
        "pressure_atm": pressure,
        "n_iterations": len(volume_chunks),
        "tau_volume_iterations": bootstrap["tau_iterations"],
        "bootstrap_block_size": bootstrap["block_size"],
        "effective_iterations": bootstrap["effective_chunks"],
        "kappa_Tmu_Pa_inverse": bootstrap["kappa"],
        "kappa_Tmu_std_Pa_inverse": bootstrap["kappa_std"],
    }

    return result, per_cycle, diagnostics


def write_csv(path, rows, fieldnames):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Saved {path}")


def validate_arguments(parser, args):
    if not 0.0 <= args.eq_frac < 1.0:
        parser.error("--eq-frac must be in [0, 1)")

    if not 0.0 < args.tail_frac <= 1.0:
        parser.error("--tail-frac must be in (0, 1]")

    if args.downsample < 1:
        parser.error("--downsample must be at least 1")

    if args.n_bootstrap < 0:
        parser.error("--n-bootstrap cannot be negative")

    low_set = args.tail_frac_low is not None
    switch_set = args.p_switch is not None

    if low_set != switch_set:
        parser.error(
            "--tail-frac-low and --p-switch must be supplied together"
        )

    if low_set and not 0.0 < args.tail_frac_low <= 1.0:
        parser.error("--tail-frac-low must be in (0, 1]")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract adsorption, volume, and compressibility data from "
            "pySIMAID output directories."
        )
    )
    parser.add_argument(
        "root_dir",
        help="Directory containing T_<temperature>/P_<pressure> directories",
    )
    parser.add_argument(
        "datafile",
        help="Empty-framework LAMMPS data file used to calculate loading",
    )
    parser.add_argument(
        "-T",
        "--temperature",
        type=float,
        required=True,
        help="Simulation temperature in K",
    )
    parser.add_argument(
        "--pmin",
        type=float,
        default=None,
        help="Minimum pressure in atm",
    )
    parser.add_argument(
        "--pmax",
        type=float,
        default=None,
        help="Maximum pressure in atm",
    )
    parser.add_argument(
        "--eq-frac",
        type=float,
        default=0.5,
        help="Fraction discarded from the start of each production run",
    )
    parser.add_argument(
        "--tail-frac",
        type=float,
        default=0.5,
        help="Fraction of final hybrid iterations included",
    )
    parser.add_argument(
        "--tail-frac-low",
        type=float,
        default=None,
        help="Tail fraction below --p-switch",
    )
    parser.add_argument(
        "--p-switch",
        type=float,
        default=None,
        help="Pressure threshold in atm for --tail-frac-low",
    )
    parser.add_argument(
        "--downsample",
        type=int,
        default=1,
        help="Retain every Nth thermo sample",
    )
    parser.add_argument(
        "--accepted-only",
        "--only-accepted-npt",
        dest="accepted_only",
        action="store_true",
        help="Include only NPT trajectories accepted by pySIMAID",
    )
    parser.add_argument(
        "--n-bootstrap",
        type=int,
        default=500,
        help="Number of block-bootstrap replicates",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=12345,
        help="Random seed used by the bootstrap",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="isotherm.csv",
        help="Per-pressure output CSV",
    )
    parser.add_argument(
        "--per-cycle-csv",
        default=None,
        help="Optional per-cycle output CSV",
    )
    parser.add_argument(
        "--bootstrap-csv",
        default=None,
        help="Optional bootstrap-diagnostics CSV",
    )

    args = parser.parse_args()
    validate_arguments(parser, args)

    pressure_directories = find_pressure_directories(
        args.root_dir,
        args.temperature,
    )

    pressure_directories = [
        (pressure, directory)
        for pressure, directory in pressure_directories
        if (args.pmin is None or pressure >= args.pmin)
        and (args.pmax is None or pressure <= args.pmax)
    ]

    if not pressure_directories:
        parser.error(
            "No matching T_<temperature>/P_<pressure> directories were found"
        )

    summaries = []
    cycles = []
    diagnostics = []

    for index, (pressure, directory) in enumerate(pressure_directories):
        processed = process_pressure(
            pressure=pressure,
            pressure_dir=directory,
            framework_datafile=args.datafile,
            temperature=args.temperature,
            equilibration_fraction=args.eq_frac,
            tail_fraction=args.tail_frac,
            tail_fraction_low=args.tail_frac_low,
            pressure_switch=args.p_switch,
            downsample=args.downsample,
            accepted_only=args.accepted_only,
            n_bootstrap=args.n_bootstrap,
            random_seed=args.seed + index,
        )

        if processed is None:
            continue

        summary, pressure_cycles, pressure_diagnostics = processed
        summaries.append(summary)
        cycles.extend(pressure_cycles)
        diagnostics.append(pressure_diagnostics)

    if not summaries:
        parser.error("No usable simulation data were extracted")

    write_csv(
        args.output,
        summaries,
        [
            "pressure_atm",
            "adsorption_mean_mmol_g",
            "adsorption_std_mmol_g",
            "volume_mean_A3",
            "volume_std_A3",
            "volume_variance_A6",
            "kappa_Tmu_Pa_inverse",
            "kappa_Tmu_std_Pa_inverse",
            "n_iterations",
            "n_adsorption_samples",
            "n_volume_samples",
        ],
    )

    if args.per_cycle_csv:
        write_csv(
            args.per_cycle_csv,
            cycles,
            [
                "pressure_atm",
                "iteration",
                "accepted",
                "N_alpha_mmol_g",
                "volume_mean_A3",
                "volume_mean_square_A6",
                "volume_variance_A6",
                "kappa_TN_Pa_inverse",
                "n_volume_samples",
            ],
        )

    if args.bootstrap_csv:
        write_csv(
            args.bootstrap_csv,
            diagnostics,
            [
                "pressure_atm",
                "n_iterations",
                "tau_volume_iterations",
                "bootstrap_block_size",
                "effective_iterations",
                "kappa_Tmu_Pa_inverse",
                "kappa_Tmu_std_Pa_inverse",
            ],
        )


if __name__ == "__main__":
    main()
