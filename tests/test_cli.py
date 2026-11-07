"""The CLI contract: exit codes, stream discipline, formats that parse.

``main`` is called in-process. A subprocess test would prove that the *installed*
entry point works — ``scripts/validate_against_scipy.py`` checks that — but it would
not give a useful traceback when a parser change breaks an argument, so this file
stays in-process and fast.
"""

from __future__ import annotations

import json

import pytest

from kinetics.cli import EXIT_NO_ANSWER, EXIT_OK, EXIT_UNREADABLE, main

MM_HEADER = ["S", "v"]
MM_ROWS = [[1, 1.0], [2, 1.8], [5, 3.4], [10, 5.0], [20, 6.6], [50, 8.0], [100, 8.7]]
INH_HEADER = ["S", "I", "v"]
INH_ROWS = [
    [2, 0, 0.73],
    [2, 4, 0.49],
    [2, 16, 0.26],
    [10, 0, 2.86],
    [10, 4, 2.11],
    [10, 16, 1.22],
    [100, 0, 8.0],
    [100, 4, 6.67],
    [100, 16, 4.44],
]
CAL_HEADER = ["concentration", "signal"]
CAL_ROWS = [[0, 0.02], [1, 0.51], [2, 1.02], [5, 2.49], [10, 4.98], [20, 10.1]]


@pytest.fixture
def mm_csv(write_csv):
    return write_csv("mm.csv", MM_HEADER, MM_ROWS)


@pytest.fixture
def inh_csv(write_csv):
    return write_csv("inhibition.csv", INH_HEADER, INH_ROWS)


@pytest.fixture
def cal_csv(write_csv):
    return write_csv("calibration.csv", CAL_HEADER, CAL_ROWS)


def test_version_flag(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    assert "kinetics 0.1.0" in capsys.readouterr().out


def test_no_subcommand_is_a_usage_error(capsys):
    with pytest.raises(SystemExit) as exc:
        main([])
    assert exc.value.code == 2


@pytest.mark.parametrize("subcommand", ["mm", "compare"])
def test_michaelis_menten_subcommands_succeed(subcommand, mm_csv, capsys):
    assert main([subcommand, mm_csv]) == EXIT_OK
    out = capsys.readouterr().out
    assert "# " in out and "vmax" in out


def test_inhibition_subcommand_succeeds(inh_csv, capsys):
    assert main(["inhibition", inh_csv]) == EXIT_OK
    out = capsys.readouterr().out
    assert "best by AICc" in out


def test_calibrate_subcommand_succeeds(cal_csv, capsys):
    assert main(["calibrate", cal_csv]) == EXIT_OK
    out = capsys.readouterr().out
    assert "LOD" in out and "LOQ" in out


@pytest.mark.parametrize("subcommand", ["mm", "compare", "inhibition", "calibrate"])
@pytest.mark.parametrize("fmt", ["markdown", "json", "html"])
def test_every_subcommand_renders_every_format(subcommand, fmt, mm_csv, inh_csv, cal_csv, capsys):
    path = {"inhibition": inh_csv, "calibrate": cal_csv}.get(subcommand, mm_csv)
    assert main([subcommand, path, "--format", fmt]) == EXIT_OK
    out = capsys.readouterr().out
    assert out.strip()
    if fmt == "json":
        json.loads(out)  # must parse, or the pipe is useless
    elif fmt == "html":
        assert out.lstrip().startswith("<!DOCTYPE html>")
    else:
        assert out.startswith("# ")


def test_json_output_of_a_fit_contains_no_nan_tokens(mm_csv, capsys):
    main(["mm", mm_csv, "--format", "json"])
    out = capsys.readouterr().out
    assert "NaN" not in out and "Infinity" not in out
    json.loads(out)


def test_out_writes_the_file_instead_of_stdout(mm_csv, tmp_path, capsys):
    target = tmp_path / "report.md"
    assert main(["mm", mm_csv, "--out", str(target)]) == EXIT_OK
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "wrote" in captured.err
    assert target.read_text(encoding="utf-8").startswith("# ")


def test_missing_file_exits_2_and_explains(tmp_path, capsys):
    assert main(["mm", str(tmp_path / "nope.csv")]) == EXIT_UNREADABLE
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "no such file" in captured.err


def test_wrong_columns_exit_2_and_name_the_columns(write_csv, capsys):
    path = write_csv("wrong.csv", ["substrate", "rate"], [[1, 2], [3, 4]])
    assert main(["inhibition", path]) == EXIT_UNREADABLE
    err = capsys.readouterr().err
    assert "missing column" in err and "I" in err


def test_non_numeric_cell_exits_2_and_names_the_line(write_csv, capsys):
    path = write_csv("bad.csv", MM_HEADER, [[1, "lots"], [2, 3]])
    assert main(["mm", path]) == EXIT_UNREADABLE
    assert "line 2" in capsys.readouterr().err


def test_header_aliases_are_accepted(write_csv, capsys):
    """Users should not have to rewrite data files to suit a parser."""
    path = write_csv(
        "aliased.csv",
        ["[S]", "Velocity"],
        [[1, 1.0], [2, 1.8], [5, 3.4], [10, 5.0], [20, 6.6], [50, 8.0]],
    )
    assert main(["mm", path]) == EXIT_OK
    assert "vmax" in capsys.readouterr().out


def test_empty_file_exits_2(write_csv, capsys):
    path = write_csv("empty.csv", MM_HEADER, [])
    assert main(["mm", path]) == EXIT_UNREADABLE


def test_readable_but_unfittable_data_exits_1(write_csv, capsys):
    """Three points cannot support a two-parameter fit with useful uncertainty."""
    path = write_csv("tiny.csv", MM_HEADER, [[1, 1.0], [2, 1.8]])
    assert main(["mm", path]) == EXIT_NO_ANSWER
    assert "no answer" in capsys.readouterr().err


def test_bootstrap_flag_records_its_seed(mm_csv, capsys):
    assert main(["mm", mm_csv, "--bootstrap", "100", "--seed", "9"]) == EXIT_OK
    out = capsys.readouterr().out
    assert "seed 9" in out and "95% CI" in out


def test_hill_flag_adds_the_second_model(mm_csv, capsys):
    assert main(["mm", mm_csv, "--hill"]) == EXIT_OK
    out = capsys.readouterr().out
    assert "hill" in out and "k_half" in out


def test_weighted_fit_uses_the_sd_column(write_csv, capsys):
    path = write_csv(
        "weighted.csv",
        ["S", "v", "sd"],
        [
            [1, 1.0, 0.1],
            [2, 1.8, 0.1],
            [5, 3.4, 0.1],
            [10, 5.0, 0.1],
            [20, 6.6, 0.1],
            [50, 8.0, 0.1],
            [100, 8.7, 0.1],
        ],
    )
    assert main(["mm", path]) == EXIT_OK
    assert "1/sd" in capsys.readouterr().out


def test_compare_ranking_excludes_the_linearisations(mm_csv, capsys):
    assert main(["compare", mm_csv]) == EXIT_OK
    out = capsys.readouterr().out
    assert "excluded from ranking" in out
    assert "different residual scale" in out


def test_a_directory_passed_as_a_dataset_exits_2(tmp_path, capsys):
    """The OSError path: exists() is true, reading is not possible."""
    assert main(["mm", str(tmp_path)]) == EXIT_UNREADABLE
    assert "cannot read" in capsys.readouterr().err


def test_a_row_with_the_wrong_number_of_fields_exits_2(write_csv, capsys):
    path = write_csv("ragged.csv", MM_HEADER, [[1, 1.0], [2, 1.8, 99], [5, 3.4]])
    assert main(["mm", path]) == EXIT_UNREADABLE
    err = capsys.readouterr().err
    assert "line 3" in err and "3 fields" in err


def test_empty_cells_are_skipped_and_reported_as_missing_columns(write_csv, capsys):
    """An sd column with a hole means the weighting was not applied, and the
    report says 'unweighted' rather than silently weighting 6 of 7 points."""
    path = write_csv(
        "hole.csv",
        ["S", "v", "sd"],
        [
            [1, 1.0, 0.1],
            [2, 1.8, 0.1],
            [5, 3.4, 0.1],
            [10, 5.0, 0.1],
            [20, 6.6, 0.1],
            [50, 8.0, 0.1],
            [100, 8.7, ""],
        ],
    )
    assert main(["mm", path]) == EXIT_OK
    assert "unweighted" in capsys.readouterr().out


def test_a_non_numeric_cell_names_the_column(write_csv, capsys):
    path = write_csv("word.csv", MM_HEADER, [[1, "trace"], [2, 1.8]])
    assert main(["mm", path]) == EXIT_UNREADABLE
    assert "v='trace'" in capsys.readouterr().err


def test_header_only_file_is_distinguished_from_an_empty_one(write_csv, capsys):
    path = write_csv("header.csv", MM_HEADER, [])
    assert main(["mm", path]) == EXIT_UNREADABLE
    assert "no measurements" in capsys.readouterr().err


def test_a_zero_byte_file_says_it_is_empty(tmp_path, capsys):
    path = tmp_path / "zero.csv"
    path.write_text("", encoding="utf-8")
    assert main(["mm", str(path)]) == EXIT_UNREADABLE
    assert "is empty" in capsys.readouterr().err
