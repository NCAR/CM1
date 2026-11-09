"""Compare CM1 domain statistics (cm1out_stats) with a stored reference.

Used as a regression check: a run whose stats differ from the reference by
more than the tolerance fails.  The reference is a plain CSV file (one row per
stats output time, one column per variable) so that a change to it shows up
in a pull request diff.

Usage:
    # compare a run with the reference (exit status 1 if they differ)
    python compare_stats.py run/cm1out_stats.ctl reference.csv

    # loose comparison (e.g. across compilers), with a floor for a near-zero variable
    python compare_stats.py run/cm1out_stats.ctl reference.csv --rtol 0.1 \
        --atol-var maxqhl=1e-5

    # write a new reference from a run
    python compare_stats.py run/cm1out_stats.ctl reference.csv --write

Only GrADS stats output (output_format = 1) is read for now.  Needs numpy only.
"""

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

# Exit statuses
PASS, DIFFER, ERROR = 0, 1, 2


def read_grads_stats(ctl_path: Path) -> dict[str, np.ndarray]:
    """Read cm1out_stats.ctl/.dat into {variable name: values over time}."""
    lines = [line.split() for line in ctl_path.read_text().splitlines() if line.strip()]
    keyword = {fields[0].lower(): fields for fields in lines}

    dset = keyword["dset"][1]
    # In GrADS, a leading "^" means the data file path is relative to the .ctl file.
    dat_path = ctl_path.parent / dset[1:] if dset.startswith("^") else Path(dset)
    options = " ".join(keyword.get("options", [])).lower()
    dtype = ">f4" if ("big_endian" in options or "byteswapped" in options) else "<f4"

    ivars = next(i for i, fields in enumerate(lines) if fields[0].lower() == "vars")
    nvars = int(lines[ivars][1])
    names = [fields[0] for fields in lines[ivars + 1 : ivars + 1 + nvars]]

    data = np.fromfile(dat_path, dtype=dtype)
    if data.size % nvars:
        raise ValueError(
            f"{dat_path}: {data.size} values is not a multiple of {nvars} variables"
        )
    ntimes = data.size // nvars
    ntimes_ctl = int(keyword["tdef"][1])
    if ntimes != ntimes_ctl:
        print(
            f"warning: {dat_path} holds {ntimes} times but {ctl_path.name} says {ntimes_ctl}"
        )

    # GrADS order: all variables at time 1, then all variables at time 2, ...
    # Keep float32 (as CM1 wrote it) so values round-trip exactly through the CSV.
    data = data.reshape(ntimes, nvars).astype(np.float32)
    return {name: data[:, i] for i, name in enumerate(names)}


def read_stats(path: Path) -> dict[str, np.ndarray]:
    if path.suffix == ".dat":
        path = path.with_suffix(".ctl")
    elif path.suffix == "":
        path = path.with_name(path.name + ".ctl")
    if path.suffix != ".ctl":
        raise ValueError(f"{path}: only GrADS stats (.ctl/.dat) are supported")
    return read_grads_stats(path)


def write_reference(stats: dict[str, np.ndarray], path: Path) -> None:
    names = list(stats)
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(names)
        # 9 significant digits round-trip a 32-bit float exactly
        for row in zip(*(stats[name] for name in names)):
            writer.writerow([f"{value:.9g}" for value in row])


def read_reference(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="") as f:
        rows = list(csv.reader(f))
    names, values = rows[0], np.array(rows[1:], dtype=np.float32)
    return {name: values[:, i] for i, name in enumerate(names)}


def first_time(mask: np.ndarray, mtime: np.ndarray) -> str:
    """Model time (s) where mask is first True, or "-" if it never is."""
    return f"{mtime[np.argmax(mask)]:.0f}" if mask.any() else "-"


def compare(
    run: dict[str, np.ndarray],
    ref: dict[str, np.ndarray],
    rtol: float,
    atol: float,
    atol_var: dict[str, float],
    flag_rtol: float,
    verbose: bool,
) -> int:
    status = PASS

    only_run = sorted(set(run) - set(ref))
    only_ref = sorted(set(ref) - set(run))
    if only_run or only_ref:
        print(f"variables only in run: {only_run}")
        print(f"variables only in reference: {only_ref}")
        status = ERROR

    ntimes_run = len(next(iter(run.values())))
    ntimes_ref = len(next(iter(ref.values())))
    if ntimes_run != ntimes_ref:
        print(f"run has {ntimes_run} output times, reference has {ntimes_ref}")
        status = ERROR
    # Keep going on a mismatch: comparing the times both have shows
    # where a run that stopped early had already gone wrong.
    nt = min(ntimes_run, ntimes_ref)
    mtime = ref.get("mtime", np.arange(nt))[:nt]

    # Columns: first time the values differ at all; first time they differ by more
    # than flag_rtol (informational only); and the largest difference, with the
    # run and reference values at that time.
    flag = f"first >{flag_rtol:.0%}"
    print(
        f"{'variable':<10} {'first diff':>10} {flag:>10} {'worst at':>9}"
        f" {'run':>12} {'reference':>12} {'abs diff':>11} {'rel diff':>9}"
    )
    ndiffer = 0
    for name in (name for name in ref if name in run):
        a = run[name][:nt].astype(np.float64)
        b = ref[name][:nt].astype(np.float64)
        absdiff = np.abs(a - b)
        # Relative to the larger magnitude, so a value near zero in either run
        # gives a relative difference of at most 1 (same sign) instead of 1e30.
        scale = np.maximum(np.abs(a), np.abs(b))
        with np.errstate(invalid="ignore", divide="ignore"):
            reldiff = np.where(scale > 0, absdiff / scale, 0.0)
        nan_mismatch = np.isnan(a) != np.isnan(b)
        differ = nan_mismatch | (absdiff > 0)
        bad = nan_mismatch | (absdiff > atol_var.get(name, atol) + rtol * scale)
        if bad.any():
            ndiffer += 1
        elif not verbose:
            continue

        if np.isnan(absdiff).all():
            iworst = 0
        else:
            iworst = int(np.nanargmax(absdiff))
        print(
            f"{name:<10} {first_time(differ, mtime):>10}"
            f" {first_time(nan_mismatch | (reldiff > flag_rtol), mtime):>10}"
            f" {mtime[iworst]:9.0f} {a[iworst]:12.4e} {b[iworst]:12.4e}"
            f" {absdiff[iworst]:11.3e} {reldiff[iworst]:9.2e}"
        )

    print(
        f"{ndiffer} of {len(ref)} variables differ beyond rtol={rtol:g}, atol={atol:g}"
        + (f" (per-variable atol: {atol_var})" if atol_var else "")
    )
    if ndiffer and status == PASS:
        status = DIFFER
    return status


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "stats", type=Path, help="cm1out_stats.ctl or .dat from the run"
    )
    parser.add_argument("reference", type=Path, help="reference CSV file")
    parser.add_argument(
        "--write",
        action="store_true",
        help="write the run's stats to the reference file",
    )
    parser.add_argument(
        "--rtol", type=float, default=0.0, help="relative tolerance (default 0: exact)"
    )
    parser.add_argument(
        "--atol", type=float, default=0.0, help="absolute tolerance (default 0: exact)"
    )
    parser.add_argument(
        "--atol-var",
        action="append",
        default=[],
        metavar="NAME=VALUE",
        help="absolute tolerance for one variable, overriding --atol (repeatable)",
    )
    parser.add_argument(
        "--flag-rtol",
        type=float,
        default=0.01,
        help="relative difference reported in the 'first >' column (default 0.01)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="list every variable, not only those that differ",
    )
    args = parser.parse_args()

    atol_var = {}
    for item in args.atol_var:
        name, _, value = item.partition("=")
        try:
            atol_var[name] = float(value)
        except ValueError:
            parser.error(f"--atol-var expects NAME=VALUE, got {item!r}")

    try:
        run = read_stats(args.stats)
        if args.write:
            write_reference(run, args.reference)
            print(
                f"wrote {len(run)} variables at {len(run['mtime'])} times to {args.reference}"
            )
            return PASS
        ref = read_reference(args.reference)
    except (OSError, ValueError, KeyError) as err:
        print(f"error: {err}")
        return ERROR

    unknown = sorted(set(atol_var) - set(ref))
    if unknown:
        print(f"error: --atol-var names not in the reference: {unknown}")
        return ERROR

    return compare(
        run, ref, args.rtol, args.atol, atol_var, args.flag_rtol, args.verbose
    )


if __name__ == "__main__":
    sys.exit(main())
