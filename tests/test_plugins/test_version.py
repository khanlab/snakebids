from __future__ import annotations

from typing import Any

import pytest
from hypothesis import given
from hypothesis import strategies as st
from pytest_mock import MockerFixture

from snakebids.plugins.version import Version, impm
from tests.helpers import allow_function_scoped


@given(version=st.text())
def test_explicit_version_used_directly(version: str):
    plug = Version(version=version)
    config: dict[Any, Any] = {}
    plug.initialize_config(config)
    assert config == {"plugins.version.version": version}


@given(version=st.text(), distribution=st.text())
@allow_function_scoped
def test_distribution_is_retrieved(
    version: str, distribution: str, mocker: MockerFixture
):
    mocker.stopall()
    mock = mocker.patch.object(impm, "version", return_value=version)
    plug = Version(distribution=distribution)
    config: dict[Any, Any] = {}
    plug.initialize_config(config)
    assert config == {"plugins.version.version": version}
    mock.assert_called_once_with(distribution)


@given(distribution=st.text())
@allow_function_scoped
def test_uninstalled_distribution_results_in_unknown(
    distribution: str, mocker: MockerFixture
):
    mocker.stopall()
    mock = mocker.patch.object(impm, "version", side_effect=impm.PackageNotFoundError)
    plug = Version(distribution=distribution)
    config: dict[Any, Any] = {}
    plug.initialize_config(config)
    assert plug.get(config, "version") == "unknown"
    mock.assert_called_once_with(distribution)


@given(version=st.text(), distribution=st.text())
def test_version_and_distribution_cannot_both_be_provided(
    version: str, distribution: str
):
    with pytest.raises(ValueError, match="version and distribution may not both"):
        Version(version=version, distribution=distribution)


def test_one_of_version_and_distribution_must_be_provided():
    with pytest.raises(ValueError, match="One of version or distribution"):
        Version()
