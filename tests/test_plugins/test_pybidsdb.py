from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Any

import pytest
from hypothesis import given
from hypothesis import strategies as st

from snakebids.plugins.pybidsdb import Pybidsdb
from snakebids.utils.utils import DEPRECATION_FLAG
from tests import strategies as sb_st


class TestAddArguments:
    def test_args_not_added_if_already_present(self):
        pybidsdb = Pybidsdb()
        parser = argparse.ArgumentParser()
        parser.add_argument("--pybidsdb-dir", dest="plugins.pybidsdb.dir")
        parser.add_argument("--pybidsdb-reset", dest="plugins.pybidsdb.reset")
        parser.add_argument("--reset-db", dest="plugins.pybidsdb.reset_db")
        pybidsdb.add_cli_arguments(parser, {})
        args = parser.parse_args(
            ["--pybidsdb-dir", "foo", "--pybidsdb-reset", "foo", "--reset-db", "foo"]
        )
        args.__dict__ = {
            "pybidsdb_dir": "foo",
            "pybidsdb_reset": "foo",
            "reset_db": "foo",
        }

    @given(title=st.text() | st.none())
    def test_args_put_in_specified_group(self, title: str | None):
        pybidsdb = Pybidsdb(title)
        parser = argparse.ArgumentParser()
        groups = {title: parser.add_argument_group(title)} if title is not None else {}
        pybidsdb.add_cli_arguments(parser, groups)

        if title is not None:
            assert groups[title]._actions
        assert parser._actions


class TestUpdateNamespace:
    @given(path=sb_st.paths())
    def test_pybidsdb_path_resolved(self, path: Path):
        pybidsdb = Pybidsdb()

        namespace = {
            f"{pybidsdb.PREFIX}.dir": path,
            f"{pybidsdb.PREFIX}.reset": False,
            f"{pybidsdb.PREFIX}.reset_db": False,
        }
        config: dict[str, Any] = {}
        pybidsdb.update_cli_namespace(namespace, config)
        # Prepare app and initial config values
        assert config["pybidsdb_dir"] == path.resolve()

    @pytest.mark.parametrize("value", [True, False])
    def test_pybidsdb_reset_and_reset_db_both_work(self, value: bool):
        pybidsdb = Pybidsdb()

        namespace = {
            f"{pybidsdb.PREFIX}.dir": Path(),
            f"{pybidsdb.PREFIX}.reset": value,
            f"{pybidsdb.PREFIX}.reset_db": not value,
        }
        config: dict[str, Any] = {}
        pybidsdb.update_cli_namespace(namespace, config)
        # Prepare app and initial config values
        assert config["pybidsdb_reset"] is True

    def test_warning_issued_for_reset_db(self, caplog: pytest.LogCaptureFixture):
        caplog.set_level(logging.WARNING)
        pybidsdb = Pybidsdb()

        namespace = {
            f"{pybidsdb.PREFIX}.dir": Path(),
            f"{pybidsdb.PREFIX}.reset": True,
            f"{pybidsdb.PREFIX}.reset_db": True,
        }
        config: dict[str, Any] = {}
        pybidsdb.update_cli_namespace(namespace, config)
        # Prepare app and initial config values
        assert "--reset-db/--reset_db will be deprecated" in caplog.text

    @pytest.mark.parametrize("value", [True, False])
    def test_namespace_cleared(self, value: bool):
        pybidsdb = Pybidsdb()

        namespace = {
            f"{pybidsdb.PREFIX}.dir": Path(),
            f"{pybidsdb.PREFIX}.reset": value,
            f"{pybidsdb.PREFIX}.reset_db": not value,
        }
        config: dict[str, Any] = {}
        pybidsdb.update_cli_namespace(namespace, config)
        # Prepare app and initial config values
        assert namespace == {}

    @pytest.mark.parametrize("value", [True, False])
    @given(path=sb_st.paths())
    def test_old_config_keys_deprecated(self, path: Path, value: bool):
        pybidsdb = Pybidsdb()

        namespace = {
            f"{pybidsdb.PREFIX}.dir": path,
            f"{pybidsdb.PREFIX}.reset": value,
            f"{pybidsdb.PREFIX}.reset_db": False,
        }
        config: dict[str, Any] = {}
        pybidsdb.update_cli_namespace(namespace, config)

        assert (
            config["pybids_db_dir"]
            == f"{DEPRECATION_FLAG}{path.resolve()}{DEPRECATION_FLAG}"
        )
        assert (
            config["pybids_db_reset"]
            == f"{DEPRECATION_FLAG}{int(value)}{DEPRECATION_FLAG}"
        )

    def test_old_key_not_set_if_pybidsdb_not_set(self):
        pybidsdb = Pybidsdb()

        namespace = {
            f"{pybidsdb.PREFIX}.dir": None,
            f"{pybidsdb.PREFIX}.reset": False,
            f"{pybidsdb.PREFIX}.reset_db": False,
        }
        config: dict[str, Any] = {}
        pybidsdb.update_cli_namespace(namespace, config)

        assert config["pybids_db_dir"] is None
