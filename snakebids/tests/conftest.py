from __future__ import annotations

import os
import tempfile
from pathlib import Path

import bids.layout
import pytest
from hypothesis import database, settings
from pyfakefs.fake_filesystem import FakeFilesystem
from upath import UPath

import snakebids.paths.resources as specs
from snakebids import resources, set_bids_spec

## Hypothesis profiles

# TODO: The default Directory Database would be much better, but it's currently
#       incompatible with pyfakefs. To fix it, we need to either patch @given or add
#       another decorator that only activates the fake filesystem within the test body

# github-actions tends to have flaky runtimes, likely due to temporary slowdowns in the
# runner, so just disable deadlines
settings.register_profile(
    "pr", deadline=None, database=database.InMemoryExampleDatabase()
)
settings.register_profile("dev", database=database.InMemoryExampleDatabase())

settings.load_profile(os.getenv("HYPOTHESIS_PROFILE", "dev"))

set_bids_spec("v0_0_0")


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
def tmpdir(request: pytest.FixtureRequest, fakefs: FakeFilesystem | None) -> Path:
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
