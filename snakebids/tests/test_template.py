from pathlib import Path

import pytest
from cookiecutter.main import cookiecutter

import snakebids
from snakebids.app import SnakeBidsApp
from snakebids.cli import SnakebidsArgs


def test_template_dry_runs_successfully(tmp_path):
    cookiecutter(
        str(Path(snakebids.__path__[0]) / "project_template"),
        no_input=True,
        output_dir=tmp_path,
    )
    app = SnakeBidsApp(
        tmp_path / "my_bids_app",
        args=SnakebidsArgs(
            force=False,
            outputdir=tmp_path / "out",
            snakemake_args=["-n"],
            args_dict={
                "bids_dir": Path("snakebids") / "tests" / "data" / "bids_bold",
                "output_dir": tmp_path / "out",
                "analysis_level": "participant",
                "smoothing_fwhm": "1",
                "filter_bold": None,
                "wildcards_bold": None,
                "path_bold": None,
                "participant_label": None,
                "exclude_participant_label": None,
            },
        ),
    )
    with pytest.raises(SystemExit) as err:
        app.run_snakemake()
    assert err.value.code == 0
