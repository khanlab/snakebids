from __future__ import annotations

import os
import tempfile
from pathlib import Path

import bids.layout
import pytest
from hypothesis import settings
from pyfakefs.fake_filesystem import FakeFilesystem

import snakebids.paths.resources as specs
from snakebids import resources

## Hypothesis profiles

# github-actions tends to have flaky runtimes, likely due to temporary slowdowns in the
# runner, so just disable deadlines
settings.register_profile("pr", deadline=None)

settings.load_profile(os.getenv("HYPOTHESIS_PROFILE", "default"))

# Fixtures


@pytest.fixture
def fakefs(
    request: pytest.FixtureRequest,
):
    """fakefs wrapper that can be disabled using the mark "disable_fakefs"

    Disabling is useful for debugging, as fakefs messes up the breakpoints in vscode

    Returns
    -------
    FakeFilesystem or None
        None if disabled
    """
    if request.node.get_closest_marker("disable_fakefs"):  # type: ignore
        return None

    return request.getfixturevalue("fs")


@pytest.fixture
def fakefs_tmpdir(
    request: pytest.FixtureRequest, fakefs: FakeFilesystem | None
) -> Path:
    """Version of tmpdir compatible with fakefs

    If fakefs is disabled, a tmpdir is returned using the builtin tmpdir fixture.
    Otherwise, an arbitrary tempdir will be created in the fakefs

    Returns
    -------
    Path
    """
    if not fakefs:
        return Path(request.getfixturevalue("tmpdir"))
    return Path(tempfile.mkdtemp())


@pytest.fixture
def bids_fs(fakefs: FakeFilesystem | None) -> FakeFilesystem | None:
    if fakefs:
        f = Path(*bids.layout.__path__, "config")
        fakefs.add_real_file(f / "bids.json")
        fakefs.add_real_file(f / "derivatives.json")
        fakefs.add_real_file(Path(*resources.__path__) / "bids_tags.json")
        fakefs.add_real_directory(Path(*specs.__path__))
    return fakefs
