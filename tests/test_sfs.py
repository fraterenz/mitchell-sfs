"""Check cell-count spectra against hand-calculated mutation matrices."""

import numpy as np
import pytest

from mitchell_sfs.sfs import SFSRaw, SiteFrequencySpectrum, site_frequency_spectrum


def test_sfs_counts_mutations_by_number_of_cells():
    raw: SFSRaw = np.array(
        [
            # Cell 1, Cell 2, Cell 3, Cell 4
            [1, 0, 0, 0],          # Mutation 1: SNV singleton
            [0, 1, 0, 0],          # Mutation 2: indel singleton
            [1, 1, 0, 0],          # Mutation 3: present in two cells
            [1, 1, 1, 1],          # Mutation 4: present in every cell
            [0, 0, 0, 0],          # Mutation 5: absent
            [0.5, 0, 1, 0],        # Mutation 6: one positive call
            [0.5, 0.5, 0.5, 0.5],  # Mutation 7: no positive calls
        ],
        dtype=np.float64,
    )
    before = raw.copy()
    expected: SiteFrequencySpectrum = {1: 3, 2: 1, 4: 1}

    result = site_frequency_spectrum(raw)
    assert result == expected
    assert all(type(k) is int and type(count) is int for k, count in result.items())

    permuted = site_frequency_spectrum(raw[::-1, ::-1])
    assert permuted == expected
    assert list(permuted) == [1, 2, 4]
    np.testing.assert_array_equal(raw, before)


@pytest.mark.parametrize(
    "raw",
    [
        np.empty((0, 4), dtype=np.float64),
        np.array([[0, 0.5], [0.5, 0]], dtype=np.float64),
    ],
    ids=["no-mutations", "no-positive-calls"],
)
def test_no_positive_calls_returns_empty_spectrum(raw):
    assert site_frequency_spectrum(raw) == {}


@pytest.mark.parametrize("call", [np.nan, np.inf, -np.inf, -1, 0.25, 2])
def test_invalid_calls_are_rejected(call):
    raw = np.array([[call, 1]], dtype=np.float64)
    with pytest.raises(ValueError, match="Mutation calls"):
        site_frequency_spectrum(raw)


@pytest.mark.parametrize(
    ("raw", "message"),
    [
        ([[1.0, 0.0]], "two-dimensional"),
        (np.array(1.0), "two-dimensional"),
        (np.ones(4, dtype=np.float64), "two-dimensional"),
        (np.ones((1, 2, 3), dtype=np.float64), "two-dimensional"),
        (np.empty((2, 0), dtype=np.float64), "sampled-cell column"),
        (np.ones((1, 2), dtype=np.int64), "dtype float64"),
        (np.ones((1, 2), dtype=np.float32), "dtype float64"),
    ],
)
def test_invalid_matrix_contract_is_rejected(raw, message):
    with pytest.raises(ValueError, match=message):
        site_frequency_spectrum(raw)
