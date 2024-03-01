from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pytest
from hypothesis import given
from pytest_mock import MockerFixture

from snakebids.app import SnakeBidsApp
from snakebids.cli import add_dynamic_args
from snakebids.tests import strategies as sb_st


class TestDeprecations:
    def test_parser_deprecated(self):
        app = SnakeBidsApp(Path(), parser=argparse.ArgumentParser())
        with pytest.warns(UserWarning, match="`SnakeBidsApp.parser` is deprecated"):
            app._check_deprecations()

    def test_config_deprecated(self):
        app = SnakeBidsApp(Path(), config={})
        with pytest.warns(UserWarning, match="`SnakeBidsApp.config` is deprecated"):
            app._check_deprecations()

    def test_args_deprecated(self):
        app = SnakeBidsApp(Path(), args=[])
        with pytest.warns(UserWarning, match="`SnakeBidsApp.args` is deprecated"):
            app._check_deprecations()

    def test_version_deprecated(self):
        app = SnakeBidsApp(Path(), version="foo")
        with pytest.warns(UserWarning, match="`SnakeBidsApp.version` is deprecated"):
            app._check_deprecations()


def test_arguments_carried_forward(mocker: MockerFixture):
    from snakebids.app import sb_plugins
    from snakebids.bidsapp import run

    mocker.stopall()
    mocker.patch.object(run, "_Runner")
    snakemake = mocker.spy(sb_plugins, "SnakemakeBidsApp")
    SnakeBidsApp(
        Path(),
        configfile_path=Path("config"),
        snakefile_path=Path("Snakefile"),
    ).run_snakemake()
    snakemake.assert_called_once_with(
        snakemake_dir=Path(),
        configfile_path=Path("config"),
        snakefile_path=Path("Snakefile"),
    )


def test_plugins_carried_forward(mocker: MockerFixture):
    from snakebids.app import bidsapp
    from snakebids.bidsapp import run

    mocker.stopall()
    mocker.patch.object(run, "_Runner")
    app_spy = mocker.spy(bidsapp, "app")
    SnakeBidsApp(
        Path(),
        configfile_path=Path("config"),
        snakefile_path=Path("Snakefile"),
        plugins=["one", "two"],  # type: ignore
    ).run_snakemake()
    app_spy.assert_called_once()
    plugins = app_spy.call_args[0][0]
    assert isinstance(plugins, list)
    assert len(plugins) == 3  # type: ignore  # noqa: PLR2004


def test_parser_can_be_directly_accessed(tmpdir: Path):
    config_path = Path(tmpdir) / "config.json"
    config_path.write_text("{}")
    app = SnakeBidsApp(
        Path(),
        configfile_path=config_path,
        snakefile_path=Path("Snakefile"),
    )
    assert app.parser is app._app.parser


def test_config_can_be_directly_accessed(tmpdir: Path):
    config_path = Path(tmpdir) / "config.json"
    config_path.write_text("{}")
    app = SnakeBidsApp(
        Path(),
        configfile_path=config_path,
        snakefile_path=Path("Snakefile"),
    )
    assert app.config is app._app.config


@given(
    parser=sb_st.everything(),
    parse_args=sb_st.everything(),
    pybids_inputs=sb_st.everything(),
)
def test_add_dynamic_args_raises_warning(
    parser: Any, parse_args: Any, pybids_inputs: Any
):
    with pytest.warns(UserWarning, match="is deprecated and no longer has any effect"):
        add_dynamic_args(parser, parse_args, pybids_inputs)


class TestGenBoutiques:
    def test_boutiques_descriptor(self, tmp_path: Path):
        configpth = tmp_path / "config.json"
        configpth.write_text("{}")
        app = SnakeBidsApp(
            Path(),
            configfile_path=configpth,
            snakefile_path=Path("Snakefile"),
            plugins=["one", "two"],  # type: ignore
        )
        descriptor_path = tmp_path / "descriptor.json"
        app.create_descriptor(descriptor_path)
        with open(descriptor_path, encoding="utf-8") as descriptor_file:
            descriptor_json = json.load(descriptor_file)
            assert "command-line" in descriptor_json
            assert "inputs" in descriptor_json
