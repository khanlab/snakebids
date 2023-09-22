from __future__ import annotations

import platform
import re
import subprocess as sp
import sys
import tempfile
from pathlib import Path
from typing import Literal, TypedDict

import copier
import more_itertools as itx
import pathvalidate
import prompt_toolkit.validation
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from typing_extensions import Unpack

if sys.version_info < (3, 11):
    import tomli as tomllib
else:
    import tomllib

import snakebids
from snakebids.tests.helpers import allow_function_scoped, needs_docker

BuildSystems = Literal["poetry", "hatch", "flit", "setuptools"]


class DataFields(TypedDict):
    full_name: str
    email: str
    app_full_name: str
    github: str
    app_description: str
    build_system: BuildSystems
    app_version: str
    create_doc_template: bool
    license: str


def get_empty_data(app_name: str, build: BuildSystems) -> DataFields:
    data: DataFields = {
        "full_name": "",
        "email": "",
        "app_full_name": app_name,
        "github": "",
        "app_description": "",
        "build_system": build,
        "app_version": "0.1.0",
        "create_doc_template": False,
        "license": "",
    }
    # poetry we need an author
    if build == "poetry":
        data["full_name"] = "John Doe"
        data["email"] = "example@email.com"
    return data


# emails complements of https://gist.github.com/cjaoude/fd9910626629b53c4d25
def invalid_emails():
    return st.sampled_from(
        [
            "plainaddress",
            "#@%^%#$@#$@#.com",
            "@example.com",
            "Joe Smith <email@example.com>",
            "email.example.com",
            "email@example@example.com",
            ".email@example.com",
            "email.@example.com",
            "email..email@example.com",
            "あいうえお@example.com",
            "email@example.com (Joe Smith)",
            "email@example",
            "email@-example.com",
            "email@example.web-",
            "email@[111.123.123.4444]",
            "email@example..com",
            "Abc..123@example.com",
            r"”(),:;<>[\]@example.com",
            "just”not”right@example.com",
            r'this\ is"really"not\allowed@example.com',
        ]
    )


@given(email=invalid_emails())
@allow_function_scoped
def test_invalid_email_raises_error(email: str, tmp_path: Path):
    data = get_empty_data("testapp", "setuptools")
    data["email"] = email
    with pytest.raises(prompt_toolkit.validation.ValidationError):
        copier.run_copy(
            str(Path(itx.first(snakebids.__path__), "project_template")),
            tmp_path / data["app_full_name"],
            data=data,
            unsafe=True,
        )


@given(
    name=st.text()
    .filter(lambda s: not re.match(r"^[a-zA-Z_][a-zA-Z_0-9]*$", s))
    .filter(lambda s: pathvalidate.is_valid_filename(s, pathvalidate.Platform.LINUX))
)
@allow_function_scoped
def test_invalid_app_name_raises_error(name: str, tmp_path: Path):
    data = get_empty_data(name, "setuptools")
    with pytest.raises(prompt_toolkit.validation.ValidationError):
        copier.run_copy(
            str(Path(itx.first(snakebids.__path__), "project_template")),
            tmp_path / data["app_full_name"],
            data=data,
            unsafe=True,
        )


@pytest.mark.parametrize(
    ["build", "build_backend"],
    [
        ("poetry", "poetry.core.masonry.api"),
        ("hatch", "hatchling.build"),
        ("flit", "flit_core.buildapi"),
        ("setuptools", "setuptools.build_meta"),
    ],
)
def test_correct_build_system_used(
    tmp_path: Path, build: BuildSystems, build_backend: str
):
    tmpdir = Path(tempfile.mkdtemp(dir=tmp_path))
    data = get_empty_data("testapp", build)
    copier.run_copy(
        str(Path(itx.first(snakebids.__path__), "project_template")),
        tmpdir / data["app_full_name"],
        data=data,
        unsafe=True,
    )
    with open(tmpdir / data["app_full_name"] / "pyproject.toml", "rb") as f:
        pyproject = tomllib.load(f)
    assert pyproject["build-system"]["build-backend"] == build_backend


@given(
    full_name=st.text(),
    email=st.emails() | st.just(""),
    app_full_name=st.from_regex(r"[a-zA-Z_][a-zA-Z_0-9]*", fullmatch=True),
    github=st.text(),
    app_description=st.text(),
    app_version=st.text(min_size=1),
    create_doc_template=st.just(False),
    license=st.text(),
)
@settings(suppress_health_check=[HealthCheck.function_scoped_fixture], deadline=1000)
@pytest.mark.parametrize(
    ["build"], [("poetry",), ("hatch",), ("flit",), ("setuptools",)]
)
def test_pyproject_correctly_formatted(
    tmp_path: Path, build: BuildSystems, **kwargs: Unpack[DataFields]
):
    tmpdir = Path(tempfile.mkdtemp(dir=tmp_path))
    kwargs["build_system"] = build
    copier.run_copy(
        str(Path(itx.first(snakebids.__path__), "project_template")),
        tmpdir / kwargs["app_full_name"],
        data=kwargs,
        unsafe=True,
    )
    with open(tmpdir / kwargs["app_full_name"] / "pyproject.toml", "rb") as f:
        pyproject = tomllib.load(f)
    if build == "poetry":
        assert kwargs["app_full_name"] == pyproject["tool"]["poetry"]["name"]
        assert kwargs["app_version"] == pyproject["tool"]["poetry"]["version"]
        assert kwargs["app_description"] == pyproject["tool"]["poetry"]["description"]
        if kwargs["license"]:
            assert kwargs["license"] == pyproject["tool"]["poetry"]["license"]
        else:
            assert "license" not in pyproject["tool"]["poetry"]
        if kwargs["full_name"]:
            assert len(pyproject["tool"]["poetry"]["authors"]) == 1
            email_tag = f" <{kwargs['email']}>" if kwargs["email"] else ""
            assert (
                f'{kwargs["full_name"]}{email_tag}'
                == pyproject["tool"]["poetry"]["authors"][0]
            )
        else:
            assert "authors" not in pyproject["tool"]["poetry"]
        return

    assert kwargs["app_full_name"] == pyproject["project"]["name"]
    assert kwargs["app_version"] == pyproject["project"]["version"]
    assert kwargs["app_description"] == pyproject["project"]["description"]
    if kwargs["license"]:
        assert kwargs["license"] == pyproject["project"]["license"]
    else:
        assert "license" not in pyproject["project"]
    if kwargs["full_name"] or kwargs["email"]:
        assert len(pyproject["project"]["authors"]) == 1
        author_obj: dict[str, str] = {}
        if kwargs["full_name"]:
            author_obj["name"] = kwargs["full_name"]
        if kwargs["email"]:
            author_obj["email"] = kwargs["email"]
        assert author_obj == pyproject["project"]["authors"][0]
    else:
        assert "authors" not in pyproject["project"]


@needs_docker(f"snakebids/test-template:{platform.python_version()}")
@pytest.mark.parametrize(
    ["build", "venv"],
    [
        ("setuptools", "setuptools"),
        ("poetry", "poetry"),
        ("hatch", "hatch"),
        ("flit", "pdm"),
    ],
)
def test_template_dry_runs_successfully(tmp_path: Path, build: BuildSystems, venv: str):
    app_name = "snakebids_app"
    data = get_empty_data(app_name, build)

    copier.run_copy(
        str(Path(itx.first(snakebids.__path__), "project_template")),
        tmp_path / app_name,
        data=data,
        unsafe=True,
    )
    cmd = sp.run(
        [
            "docker",
            "run",
            "-v",
            f"{tmp_path / app_name}:/app",
            "--rm",
            f"snakebids/test-template:{platform.python_version()}",
            venv,
            app_name,
        ],
        capture_output=True,
        check=False,
    )
    try:
        cmd.check_returncode()
    except Exception as err:
        print(cmd.stdout)
        print(cmd.stderr, file=sys.stderr)
        raise err
    assert "All set" in cmd.stdout.decode()
