from __future__ import annotations

from io import StringIO
from pathlib import Path

import pytest
from hypothesis import given
from hypothesis import strategies as st
from pyfakefs.fake_filesystem import FakeFilesystem
from pytest_mock import MockerFixture
from pytest_mock.plugin import MockType

import snakebids.io.config as configio
import snakebids.io.yaml as yamlio
from snakebids.tests.helpers import allow_function_scoped


@given(path=st.text(st.characters(blacklist_characters=["\x85"]), min_size=1).map(Path))
def test_paths_formatted_as_str(path: Path):
    string = StringIO()
    yaml = yamlio.get_yaml_io()
    yaml.dump({"key": path}, string)
    string.seek(0, 0)
    assert yaml.load(string)["key"] == str(path)


class TestWriteConfig:
    def io_mocks(self, mocker: MockerFixture) -> dict[str, MockType]:
        return {
            "mopen": mocker.patch.object(configio, "open", mocker.mock_open()),
            "jsondump": mocker.patch.object(configio.json, "dump"),
            "mkdir": mocker.patch.object(configio.Path, "mkdir"),
            "yamldump": mocker.patch.object(yamlio.YAML, "dump"),
        }

    @allow_function_scoped
    @given(
        ext=st.sampled_from([".json", ".yaml", ".yml"]),
        path=st.text().map(Path).filter(lambda p: p != Path() and not p.exists()),
    )
    def test_writes_correct_format(self, ext: str, path: Path, mocker: MockerFixture):
        mocker.stopall()
        mocks = self.io_mocks(mocker)
        path = path.with_suffix(ext)
        configio.write_config(path, {})
        if ext == ".json":
            mocks["mopen"].assert_called_once_with(path, "w", encoding="utf-8")
            mocks["jsondump"].assert_called_once()
            mocks["yamldump"].assert_not_called()
        else:
            mocks["mopen"].assert_not_called()
            mocks["jsondump"].assert_not_called()
            mocks["yamldump"].assert_called_once_with({}, path)

    @allow_function_scoped
    @given(path=st.text().filter(lambda s: s not in {".", ""}))
    def test_doesnt_overwrite_file(self, path: str, fakefs: FakeFilesystem):
        fakefs.reset()
        fakefs.create_file(path)
        with pytest.raises(FileExistsError, match="already exists"):
            configio.write_config(path, {})

    @allow_function_scoped
    @given(path=st.text().filter(lambda s: s not in {".", ""}))
    def test_overwrites_file_if_forced(
        self, path: str, fakefs: FakeFilesystem, mocker: MockerFixture
    ):
        mocker.stopall()
        fakefs.reset()
        fakefs.create_file(path)
        mocks = self.io_mocks(mocker)
        configio.write_config(path, {}, force_overwrite=True)
        assert mocks["jsondump"].call_count ^ mocks["yamldump"].call_count
