"""Export per-donor site-frequency spectra with donor ages and sample sizes."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from mitchell_sfs.sfs import SFSRaw, SiteFrequencySpectrum, site_frequency_spectrum


def donor_metadata() -> pd.DataFrame:
    """Look up donor ages included in per-donor exports.

    Returns
    -------
    pd.DataFrame
        Built-in donor records with columns ``name``, ``age`` (integer
        years), ``cells``, and ``clones``. The latter two are reference
        metadata; the SFS sample size comes from the input matrix.
    """
    return pd.DataFrame.from_records(
        [
            {"name": "CB001", "age": 0, "cells": 216, "clones": 0},
            {"name": "CB002", "age": 0, "cells": 390, "clones": 0},
            {"name": "KX001", "age": 29, "cells": 407, "clones": 0},
            {"name": "KX002", "age": 38, "cells": 380, "clones": 1},
            {"name": "SX001", "age": 48, "cells": 362, "clones": 0},
            {"name": "AX001", "age": 63, "cells": 361, "clones": 1},
            {"name": "KX008", "age": 76, "cells": 367, "clones": 12},
            {"name": "KX004", "age": 77, "cells": 451, "clones": 15},
            {"name": "KX003", "age": 81, "cells": 328, "clones": 13},
        ]
    )


def load_sfs_raw(path: Path) -> SFSRaw:
    """Read mutation calls with the same mutation-by-cell orientation as the CSV.

    Parameters
    ----------
    path
        CSV whose first column contains mutation identifiers and whose
        remaining columns contain sampled-cell calls.

    Returns
    -------
    SFSRaw
        Float64 array with mutation rows and cell columns. Identifier and
        column labels are omitted. The numerical function validates the
        shape and call values before computing the spectrum.

    Raises
    ------
    OSError
        If the input cannot be read.
    ValueError
        If the CSV cannot be parsed or its calls cannot be converted to
        floating-point values.
    """
    matrix = pd.read_csv(path, index_col=0)
    return matrix.to_numpy(dtype=np.float64)


def write_sfs_csv(
    spectrum: SiteFrequencySpectrum,
    path: Path,
    *,
    donor: str,
    age_years: int,
    sampled_cells: int,
) -> None:
    """Write all positive cell-count bins, filling unobserved bins with zeros.

    Parameters
    ----------
    spectrum
        Sparse spectrum computed from a matrix with ``sampled_cells``
        columns. Its keys must lie between 1 and ``sampled_cells``.
    path
        Destination CSV in an existing directory. Existing contents are
        replaced.
    donor
        Donor identifier written on every output row.
    age_years
        Donor age in integer years, written on every output row.
    sampled_cells
        Actual number of cell columns in the input matrix, at least one.
        The caller supplies a spectrum aligned with this sample size.

    Raises
    ------
    OSError
        If the output cannot be written.

    Notes
    -----
    Columns are ``donor``, ``age_years``, ``sampled_cells``,
    ``cells_with_mutation``, and ``mutation_count``. Rows cover bins 1
    through ``sampled_cells``. An empty spectrum yields zero mutation
    counts in all those bins.
    """
    bins = range(1, sampled_cells + 1)
    pd.DataFrame(
        {
            "donor": donor,
            "age_years": age_years,
            "sampled_cells": sampled_cells,
            "cells_with_mutation": list(bins),
            "mutation_count": [spectrum.get(k, 0) for k in bins],
        }
    ).to_csv(path, index=False)


def main(argv: list[str] | None = None) -> None:
    """Export spectra for discovered donors using positional input and output paths.

    Parameters
    ----------
    argv
        Command-line arguments excluding the executable name. ``None``
        reads the process arguments. Input defaults to
        ``./Mitchell et al 2022 - Sanger`` and output to ``INPUTDIR/sfs``.

    Raises
    ------
    SystemExit
        With status 2 for invalid arguments, missing donor metadata, or
        input/output failures; with status 0 when displaying help.
    """
    parser = argparse.ArgumentParser(
        description="Compute donor SFS using both SNVs and indels."
    )
    parser.add_argument(
        "input_dir",
        metavar="INPUTDIR",
        nargs="?",
        type=Path,
        default=Path("Mitchell et al 2022 - Sanger"),
        help="Input directory (default: %(default)s).",
    )
    parser.add_argument(
        "output_dir",
        metavar="OUTPUTDIR",
        nargs="?",
        type=Path,
        help="Output directory (default: INPUTDIR/sfs).",
    )
    args = parser.parse_args(argv)

    if not args.input_dir.is_dir():
        parser.error(
            "Input directory does not exist or is not a directory: "
            f"{args.input_dir}"
        )

    inputs = {
        path.stem.removeprefix("mutMatrix"): path
        for path in sorted(args.input_dir.glob("mutMatrix*.csv"))
        if path.is_file()
    }
    if not inputs:
        parser.error(f"No mutMatrix*.csv files found in {args.input_dir}")

    metadata = donor_metadata().set_index("name")
    unknown = sorted(set(inputs) - set(metadata.index))
    if unknown:
        parser.error(f"Missing donor metadata: {', '.join(unknown)}")

    output_dir = args.output_dir or args.input_dir / "sfs"
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
    except OSError as error:
        parser.error(f"Cannot create output directory {output_dir}: {error}")

    for donor, input_path in inputs.items():
        age_years = int(metadata.at[donor, "age"])
        output_path = output_dir / f"{donor}.csv"

        try:
            raw = load_sfs_raw(input_path)
            spectrum = site_frequency_spectrum(raw)
            write_sfs_csv(
                spectrum,
                output_path,
                donor=donor,
                age_years=age_years,
                sampled_cells=raw.shape[1],
            )
        except (OSError, ValueError) as error:
            parser.error(f"{donor}: {error}")

        print(f"{donor}: wrote {output_path}")


if __name__ == "__main__":
    main()
