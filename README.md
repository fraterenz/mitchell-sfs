# mitchell-sfs

Compute a site-frequency spectrum (SFS) for each donor in the Mitchell et al.
(2022) mutation matrices: the number of mutations present in exactly `K`
sampled cells. All mutation types, including SNVs and indels, are included;
`mutType<DONOR>.csv` annotations are not required.

## Setup and usage

Install [uv](https://docs.astral.sh/uv/getting-started/installation/), then run
from this project directory (Python 3.14 is pinned):

```bash
uv run mitchell-sfs
uv run mitchell-sfs "/path/to/input" "/path/to/results"
uv run mitchell-sfs --help
```

Both positional arguments are optional: `INPUTDIR` defaults to
`./Mitchell et al 2022 - Sanger` relative to the working directory, and
`OUTPUTDIR` to `<INPUTDIR>/sfs`. Quote paths containing spaces.

## Inputs and interpretation

The command processes `mutMatrix<DONOR>.csv` files (e.g. `mutMatrixCB001.csv`)
directly in the input directory, one donor at a time.

- Each row is a mutation; the first column holds its identifier, and the
  remaining columns represent sampled cells.
- Calls must be `0`, `0.5`, or `1`. A call of `1` means present; `0` and `0.5`
  are treated as absent for this calculation.
- Missing or nonfinite calls are rejected.

The SFS uses sampled-cell counts, without frequency normalization or folding.

Donor ages come from the metadata table in `src/mitchell_sfs/main.py`;
unrecognized donor IDs cause an error before export. Sample size comes from
the matrix, not the metadata table's reference `cells` field.

## Outputs

Each donor produces `<DONOR>.csv` (e.g. `CB001.csv`), replacing any existing
file for that donor.

| Column | Meaning |
| --- | --- |
| `donor` | Donor ID extracted from the input filename |
| `age_years` | Donor age in integer years, from the built-in metadata table |
| `sampled_cells` | Total number of cell columns in the input matrix |
| `cells_with_mutation` | Number of sampled cells containing a mutation (`K`) |
| `mutation_count` | Number of mutations present in exactly `K` cells |

Output includes bins `1` through `sampled_cells`, including mutations present
in every cell and zero counts for unobserved bins. Mutations absent from all
cells (bin `0`) are excluded.

A **synthetic four-cell example**, using the CB001 donor label, would be:

```csv
donor,age_years,sampled_cells,cells_with_mutation,mutation_count
CB001,0,4,1,3
CB001,0,4,2,1
CB001,0,4,3,0
CB001,0,4,4,1
```

## SFS plots

Generate two PDFs from the exported per-donor CSVs:

```bash
uv run python -m mitchell_sfs.plot
uv run python -m mitchell_sfs.plot "/path/to/sfs"
```

Both PDFs are saved in the SFS directory: `sfs_by_donor.pdf` and
`sfs_by_donor_normalised.pdf`. Donor panels are ordered by age, use log–log
scales, share x and y axes within each figure, and omit zero-count bins.
The normalised figure divides cell counts by each donor's `sampled_cells`
and mutation counts by that donor's total mutation count (so they sum to 1).
Both figures use the same plotting function, with a `normalise` option.

## Python interface

```python
import numpy as np

from mitchell_sfs.sfs import SFSRaw, SiteFrequencySpectrum, site_frequency_spectrum

# Rows are mutations; columns are sampled cells.
raw: SFSRaw = np.array(
    [[1, 0], [1, 1], [0.5, 0]],
    dtype=np.float64,
)
spectrum: SiteFrequencySpectrum = site_frequency_spectrum(raw)
assert spectrum == {1: 1, 2: 1}
```

`SFSRaw` is `Float64[np.ndarray, "mutation cell"]`; `SiteFrequencySpectrum`
is `dict[K, SFS]`, where both aliases are `int`. Unlike CSV output, the function
returns only observed positive bins (`{}` for an empty mutation set). It checks
shape, dtype, and calls at runtime and leaves the input unchanged.

## Tests

```bash
uv run pytest tests/test_sfs.py tests/test_cli.py
```

Tests use synthetic matrices and temporary CSVs to cover calculation,
validation, donor metadata, and output bins.

## Code and data licensing

The code is licensed under MIT; see [LICENSE](LICENSE).

The dataset requires separate authorized access and is not covered by the
software license. The provider's data-use terms govern its use and the
redistribution of data and derived outputs.

Git allowlists project metadata, Python source, and synthetic tests; data,
results, environments, and archives are ignored, even inside the project.
