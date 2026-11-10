"""Plot raw and per-donor normalised SFS panels from exported CSVs."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def plot_spectra(
    spectra: list[pd.DataFrame], output: Path, *, normalise: bool = False
) -> None:
    """Save shared-axis, log–log donor panels, omitting zero-count bins."""
    columns = min(3, len(spectra))
    rows = (len(spectra) + columns - 1) // columns
    fig, axes = plt.subplots(
        rows, columns, figsize=(4 * columns, 3.2 * rows),
        sharex=True, sharey=True, squeeze=False, layout="constrained",
    )
    for ax, spectrum in zip(axes.flat, spectra):
        metadata = spectrum.iloc[0]
        x = spectrum["cells_with_mutation"].to_numpy(dtype=float)
        y = spectrum["mutation_count"].to_numpy(dtype=float)
        positive = (x > 0) & (y > 0)
        if normalise:
            x = x / metadata["sampled_cells"]
            y = y / y.sum()
        ax.scatter(x[positive], y[positive], s=9, color="#2563a6", alpha=0.8,
                   linewidths=0)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(
            f"{metadata['donor']} · age {metadata['age_years']} years\n"
            f"{metadata['sampled_cells']} cells", fontsize=11,
        )
        ax.grid(which="major", alpha=0.25, linewidth=0.6)
        ax.set_axisbelow(True)
    for ax in list(axes.flat)[len(spectra):]:
        ax.set_visible(False)
    fig.supxlabel(
        "Fraction of sequenced cells carrying mutation (k / N)"
        if normalise else "Number of cells carrying mutation (k)"
    )
    fig.supylabel(
        "Normalised SFS (fraction of mutations)"
        if normalise else "SFS (number of mutations)"
    )
    fig.suptitle("Mitchell et al. (2022) — " + (
        "normalised site-frequency spectra" if normalise else "site-frequency spectra"
    ))
    fig.savefig(output)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "sfsdir", type=Path, nargs="?",
        default=Path("Mitchell et al 2022 - Sanger/sfs"),
        help="Directory containing per-donor SFS CSVs; PDFs are saved here.",
    )
    args = parser.parse_args()
    spectra = [pd.read_csv(path) for path in sorted(args.sfsdir.glob("*.csv"))]
    if not spectra:
        parser.error(f"No SFS CSV files found in {args.sfsdir}")
    spectra.sort(key=lambda data: (data.iloc[0]["age_years"], data.iloc[0]["donor"]))
    for normalise, filename in (
        (False, "sfs_by_donor.pdf"),
        (True, "sfs_by_donor_normalised.pdf"),
    ):
        output = args.sfsdir / filename
        plot_spectra(spectra, output, normalise=normalise)
        print(output)


if __name__ == "__main__":
    main()
