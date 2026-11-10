"""Count mutations by the number of sampled cells in which they are present."""

import numpy as np
from jaxtyping import Float64


type SFSRaw = Float64[np.ndarray, "mutation cell"]

# K is a sampled-cell count; SFS is the mutation count in that bin.
type K = int
type SFS = int
type SiteFrequencySpectrum = dict[K, SFS]


def site_frequency_spectrum(raw: SFSRaw) -> SiteFrequencySpectrum:
    """Count mutations present in exactly each positive number of sampled cells.

    Parameters
    ----------
    raw
        Two-dimensional float64 array with mutation rows and sampled-cell
        columns. Calls must be 0, 0.5, or 1. Only 1 counts as present;
        0 and 0.5 count as absent. All mutation types, including SNVs and
        indels, are considered. The input is not modified.

    Returns
    -------
    SiteFrequencySpectrum
        Dictionary mapping positive cell counts to mutation counts, with
        ascending keys and ordinary Python integer keys and values.
        Unobserved bins and the zero-cell bin are omitted. Mutations present
        in every sampled cell are retained. An empty mutation set, or one
        without positive calls, returns an empty dictionary.

    Raises
    ------
    ValueError
        If the input is not a two-dimensional float64 NumPy array, has no
        cell columns, or contains calls outside 0, 0.5, and 1, including
        missing or nonfinite values.

    Notes
    -----
    These are sampled-cell counts, without frequency normalization or
    folding. The shape and dtype annotation is supplemented by explicit
    runtime checks.
    """
    if not isinstance(raw, np.ndarray) or raw.ndim != 2:
        raise ValueError(
            "raw must be a two-dimensional NumPy array with mutation rows "
            "and cell columns."
        )
    if raw.dtype != np.dtype(np.float64):
        raise ValueError("raw must have dtype float64.")
    if raw.shape[1] == 0:
        raise ValueError("raw must have at least one sampled-cell column.")
    if not np.isin(raw, (0.0, 0.5, 1.0)).all():
        raise ValueError(
            "Mutation calls must be 0, 0.5, or 1; missing and nonfinite "
            "values are not allowed."
        )

    cells_with_mutation = np.count_nonzero(raw == 1, axis=1)
    cell_counts, mutation_counts = np.unique(
        cells_with_mutation[cells_with_mutation > 0],
        return_counts=True,
    )
    return {
        int(k): int(count)
        for k, count in zip(cell_counts, mutation_counts, strict=True)
    }
