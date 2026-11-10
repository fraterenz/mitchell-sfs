"""Check donor discovery and complete CSV exports using synthetic inputs."""

from pathlib import Path

import pandas as pd
import pytest

from mitchell_sfs.main import main


def _write_matrix(path: Path, rows: list[list[float]], *, cells: int) -> None:
    pd.DataFrame(
        rows,
        index=[f"mutation{i + 1}" for i in range(len(rows))],
        columns=[f"cell{i + 1}" for i in range(cells)],
        dtype=float,
    ).to_csv(path)


@pytest.fixture
def input_dir(tmp_path):
    path = tmp_path / "Mitchell et al 2022 - Sanger"
    path.mkdir()
    _write_matrix(
        path / "mutMatrixCB001.csv",
        [[1, 0, 0, 0], [1, 1, 0, 0], [0.5, 0, 0, 0]],
        cells=4,
    )
    pd.DataFrame(
        {"x": ["SNV", "INDEL", "SNV"]}, index=[1, 2, 3]
    ).to_csv(path / "mutTypeCB001.csv")
    _write_matrix(
        path / "mutMatrixKX001.csv",
        [[1, 1, 1], [1, 0.5, 0]],
        cells=3,
    )
    return path


@pytest.mark.parametrize("invocation", ["defaults", "input-only", "explicit-output"])
def test_cli_exports_complete_bins_with_donor_ages(
    input_dir, tmp_path, monkeypatch, invocation
):
    monkeypatch.chdir(tmp_path)
    output_dir = input_dir / "sfs"
    if invocation == "defaults":
        argv = []
    elif invocation == "input-only":
        argv = [str(input_dir)]
    else:
        output_dir = tmp_path / "export folder"
        argv = [str(input_dir), str(output_dir)]

    main(argv)

    assert sorted(path.name for path in output_dir.iterdir()) == [
        "CB001.csv",
        "KX001.csv",
    ]
    pd.testing.assert_frame_equal(
        pd.read_csv(output_dir / "CB001.csv"),
        pd.DataFrame(
            {
                "donor": ["CB001"] * 4,
                "age_years": [0] * 4,
                "sampled_cells": [4] * 4,
                "cells_with_mutation": [1, 2, 3, 4],
                "mutation_count": [1, 1, 0, 0],
            }
        ),
    )
    pd.testing.assert_frame_equal(
        pd.read_csv(output_dir / "KX001.csv"),
        pd.DataFrame(
            {
                "donor": ["KX001"] * 3,
                "age_years": [29] * 3,
                "sampled_cells": [3] * 3,
                "cells_with_mutation": [1, 2, 3],
                "mutation_count": [1, 0, 1],
            }
        ),
    )


@pytest.mark.parametrize("rows", [[], [[0, 0.5, 0, 0.5]]])
def test_empty_spectrum_retains_sample_metadata(tmp_path, rows):
    _write_matrix(tmp_path / "mutMatrixKX001.csv", rows, cells=4)

    main([str(tmp_path)])

    pd.testing.assert_frame_equal(
        pd.read_csv(tmp_path / "sfs" / "KX001.csv"),
        pd.DataFrame(
            {
                "donor": ["KX001"] * 4,
                "age_years": [29] * 4,
                "sampled_cells": [4] * 4,
                "cells_with_mutation": [1, 2, 3, 4],
                "mutation_count": [0] * 4,
            }
        ),
    )


def test_unknown_donor_is_rejected_before_output_creation(tmp_path, capsys):
    _write_matrix(tmp_path / "mutMatrixUNKNOWN.csv", [[1]], cells=1)
    output_dir = tmp_path / "sfs"

    with pytest.raises(SystemExit) as error:
        main([str(tmp_path), str(output_dir)])

    assert error.value.code == 2
    assert "Missing donor metadata: UNKNOWN" in capsys.readouterr().err
    assert not output_dir.exists()


@pytest.mark.parametrize("exists", [False, True], ids=["missing", "empty"])
def test_missing_or_empty_input_directory_is_reported(tmp_path, capsys, exists):
    input_dir = tmp_path / "input"
    if exists:
        input_dir.mkdir()

    with pytest.raises(SystemExit) as error:
        main([str(input_dir)])

    assert error.value.code == 2
    message = capsys.readouterr().err
    expected = "No mutMatrix*.csv" if exists else "Input directory does not exist"
    assert expected in message


def test_invalid_calls_are_reported_with_donor_context(tmp_path, capsys):
    _write_matrix(tmp_path / "mutMatrixCB001.csv", [[2, 0]], cells=2)

    with pytest.raises(SystemExit) as error:
        main([str(tmp_path)])

    assert error.value.code == 2
    assert "CB001: Mutation calls" in capsys.readouterr().err
    assert not (tmp_path / "sfs" / "CB001.csv").exists()


def test_output_directory_failure_is_reported(input_dir, tmp_path, capsys):
    output_path = tmp_path / "occupied"
    output_path.write_text("This is a file, not a directory.")

    with pytest.raises(SystemExit) as error:
        main([str(input_dir), str(output_path)])

    assert error.value.code == 2
    assert "Cannot create output directory" in capsys.readouterr().err
