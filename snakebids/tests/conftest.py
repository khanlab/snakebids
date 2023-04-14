# pylint: disable=redefined-outer-name
import tempfile
from pathlib import Path
from typing import Optional

import bids.layout
import pytest
from pyfakefs.fake_filesystem import FakeFilesystem


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
    if request.node.get_closest_marker("disable_fakefs"):
        return None

    return request.getfixturevalue("fs")


@pytest.fixture
def fakefs_tmpdir(
    request: pytest.FixtureRequest, fakefs: Optional[FakeFilesystem]
) -> Path:
    """Version of tmpdir compatible with fakefs

    If fakefs is disabled, a tmpdir is returned using the builtin tmpdir fixture.
    Otherwise, an arbitrary tempdir will be created in the fakefs

    Returns
    -------
    Path
    """
    if not fakefs:
        return request.getfixturevalue("tmpdir")
    return Path(tempfile.mkdtemp())


@pytest.fixture
def bids_fs(fakefs: Optional[FakeFilesystem]) -> FakeFilesystem | None:
    if fakefs:
        f = Path(*bids.layout.__path__, "config")
        fakefs.add_real_file(f / "bids.json")
        fakefs.add_real_file(f / "derivatives.json")
    return fakefs
