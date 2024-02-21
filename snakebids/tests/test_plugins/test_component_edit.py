from __future__ import annotations

import argparse
import copy
import functools as ft
import itertools as it
import re
from pathlib import Path
from typing import Any, ClassVar

import pytest
from hypothesis import given
from hypothesis import strategies as st

import snakebids.tests.strategies as sb_st
from snakebids.exceptions import MisspecifiedCliFilterError
from snakebids.plugins.component_edit import ComponentEdit
from snakebids.tests.helpers import allow_function_scoped
from snakebids.types import InputsConfig, OptionalFilter


class TestAddArguments:
    mock_args_special: ClassVar[list[str]] = ["--derivatives", "path/to/nowhere"]
    mock_basic_args: ClassVar[list[str]] = [
        "script_name",
        "path/to/input",
        "path/to/output",
        "participant",
    ]
    mock_all_args: ClassVar[list[str]] = mock_basic_args + mock_args_special

    @given(key=st.text())
    def test_works_when_key_not_present(self, key: str):
        comp_edit = ComponentEdit(components_key=key)
        p = argparse.ArgumentParser()
        comp_edit.add_cli_arguments(p, {})
        assert len(p._actions) == 1

    @given(key=st.text())
    @allow_function_scoped
    def test_config_key_can_be_set(self, key: str):
        comp_edit = ComponentEdit(components_key=key)
        p = argparse.ArgumentParser()
        comp_edit.add_cli_arguments(p, {key: {"comp": {}}})
        assert len(p._actions) == 4  # noqa: PLR2004

    @given(sb_st.inputs_configs())
    @allow_function_scoped
    def test_dynamic_inputs(self, pybids_inputs: InputsConfig):
        comp_edit = ComponentEdit()
        p = argparse.ArgumentParser()
        comp_edit.add_cli_arguments(p, {"pybids_inputs": pybids_inputs})
        magic_filters = list(
            it.chain.from_iterable(
                [[f"--filter-{key}", "entity=value"] for key in pybids_inputs]
            )
        )
        magic_wildcards = list(
            it.chain.from_iterable(
                [[f"--wildcards-{key}", "test"] for key in pybids_inputs]
            )
        )
        magic_path = list(
            it.chain.from_iterable([[f"--path-{key}", "test"] for key in pybids_inputs])
        )
        argv = magic_filters + magic_wildcards + magic_path
        args = p.parse_args(argv)
        for key in pybids_inputs:
            assert args.__dict__[f"{comp_edit.PREFIX}.path.{key}"] == "test"
            assert args.__dict__[f"{comp_edit.PREFIX}.filter.{key}"] == {
                "entity": "value"
            }
            assert args.__dict__[f"{comp_edit.PREFIX}.wildcards.{key}"] == ["test"]

    @given(
        pybids_inputs=sb_st.inputs_configs(),
        flag=st.from_regex(
            re.compile(r"(?:required)|(?:any)", re.IGNORECASE), fullmatch=True
        ),
    )
    @allow_function_scoped
    def test_required_filters(self, pybids_inputs: InputsConfig, flag: str):
        p = argparse.ArgumentParser()
        comp_edit = ComponentEdit()
        comp_edit.add_cli_arguments(p, {"pybids_inputs": pybids_inputs})
        argv = list(
            it.chain.from_iterable(
                [[f"--filter-{key}", f"entity:{flag}"] for key in pybids_inputs]
            )
        )

        args = p.parse_args(argv)
        for key in pybids_inputs:
            assert args.__dict__[f"{comp_edit.PREFIX}.filter.{key}"]["entity"] is True

    @given(
        pybids_inputs=sb_st.inputs_configs(),
        flag=st.from_regex(re.compile(r"optional", re.IGNORECASE), fullmatch=True),
    )
    @allow_function_scoped
    def test_optional_filters(self, pybids_inputs: InputsConfig, flag: str):
        p = argparse.ArgumentParser()
        comp_edit = ComponentEdit()
        comp_edit.add_cli_arguments(p, {"pybids_inputs": pybids_inputs})
        argv = list(
            it.chain.from_iterable(
                [[f"--filter-{key}", f"entity:{flag}"] for key in pybids_inputs]
            )
        )

        args = p.parse_args(argv)
        for key in pybids_inputs:
            assert (
                args.__dict__[f"{comp_edit.PREFIX}.filter.{key}"]["entity"]
                is OptionalFilter
            )

    @given(
        pybids_inputs=sb_st.inputs_configs(),
        flag=st.from_regex(re.compile(r"none", re.IGNORECASE), fullmatch=True),
    )
    @allow_function_scoped
    def test_none_filters(self, pybids_inputs: InputsConfig, flag: str):
        p = argparse.ArgumentParser()
        comp_edit = ComponentEdit()
        comp_edit.add_cli_arguments(p, {"pybids_inputs": pybids_inputs})
        argv = list(
            it.chain.from_iterable(
                [[f"--filter-{key}", f"entity:{flag}"] for key in pybids_inputs]
            )
        )

        args = p.parse_args(argv)
        for key in pybids_inputs:
            assert args.__dict__[f"{comp_edit.PREFIX}.filter.{key}"]["entity"] is False

    @given(
        pybids_inputs=sb_st.inputs_configs(min_size=1, max_size=1),
        flag=st.text().filter(
            lambda s: s.lower() not in {"none", "any", "optional", "required"}
        ),
    )
    def test_filter_with_bad_flag_errors(self, pybids_inputs: InputsConfig, flag: str):
        p = argparse.ArgumentParser()
        comp_edit = ComponentEdit()
        comp_edit.add_cli_arguments(p, {"pybids_inputs": pybids_inputs})
        argv = list(
            it.chain.from_iterable(
                [[f"--filter-{key}", f"entity:{flag}"] for key in pybids_inputs]
            )
        )

        with pytest.raises(
            MisspecifiedCliFilterError,
            match=re.compile(rf"following filter provided.*entity:{re.escape(flag)}"),
        ):
            p.parse_args(argv)

    @given(
        pybids_inputs=sb_st.inputs_configs(min_size=1, max_size=1),
        filt=st.text().filter(
            lambda s: "=" not in s and ":" not in s and not s.startswith("-")
        ),
    )
    def test_filters_with_no_delimiter_errors(
        self, pybids_inputs: InputsConfig, filt: str
    ):
        p = argparse.ArgumentParser()
        comp_edit = ComponentEdit()
        comp_edit.add_cli_arguments(p, {"pybids_inputs": pybids_inputs})
        argv = list(
            it.chain.from_iterable([[f"--filter-{key}", filt] for key in pybids_inputs])
        )

        with pytest.raises(
            MisspecifiedCliFilterError,
            match=re.compile(rf"following filter provided.*{re.escape(filt)}"),
        ):
            p.parse_args(argv)


@st.composite
def two_configs(
    draw: st.DrawFn, filters: bool = True, wildcards: bool = True, paths: bool = True
) -> tuple[InputsConfig, InputsConfig]:
    inputs_configs = ft.partial(
        sb_st.inputs_configs, filters=filters, wildcards=wildcards, paths=paths
    )
    first = draw(inputs_configs())
    if not first:
        return first, {}
    second = draw(inputs_configs(keys=st.sampled_from(list(first.keys()))))
    return first, second


class TestUpdateNamespace:
    def make_namespace(self, config: InputsConfig):
        namespace: dict[str, Any] = {}
        for comp, conf in config.items():
            if "filters" in conf:
                namespace[f"{ComponentEdit.PREFIX}.filter.{comp}"] = conf["filters"]
            if "wildcards" in conf:
                namespace[f"{ComponentEdit.PREFIX}.wildcards.{comp}"] = conf[
                    "wildcards"
                ]
            if "custom_path" in conf:
                namespace[f"{ComponentEdit.PREFIX}.path.{comp}"] = conf["custom_path"]
        return namespace

    def test_works_when_key_not_present(self):
        comp_edit = ComponentEdit()
        config: dict[Any, Any] = {}
        comp_edit.update_cli_namespace({}, config)
        assert config == {}

    @given(configs=two_configs(wildcards=False, paths=False))
    def test_filters_are_joined_with_preference_to_namespace(
        self, configs: tuple[InputsConfig, InputsConfig]
    ):
        comp_edit = ComponentEdit()
        config = configs[0]
        namespace_orig = configs[1]
        namespace = self.make_namespace(configs[1])
        comp_edit.update_cli_namespace(namespace, {"pybids_inputs": config})
        for comp, conf in namespace_orig.items():
            for entity, filter_ in conf.get("filters", {}).items():
                assert config[comp].get("filters", {})[entity] == filter_

    @given(configs=two_configs(wildcards=False, paths=False))
    def test_original_filters_remain_intact(
        self, configs: tuple[InputsConfig, InputsConfig]
    ):
        comp_edit = ComponentEdit()
        config = configs[0]
        config_orig = copy.deepcopy(config)
        namespace_orig = configs[1]
        namespace = self.make_namespace(configs[1])
        comp_edit.update_cli_namespace(namespace, {"pybids_inputs": config})
        for comp, conf in config_orig.items():
            for entity, filter_ in conf.get("filters", {}).items():
                if entity not in namespace_orig.get(comp, {}).get("filters", {}):
                    assert config[comp].get("filters", {})[entity] == filter_

    @given(config=sb_st.inputs_configs(wildcards=False, paths=False))
    def test_optional_filter_removes_filters(self, config: InputsConfig):
        comp_edit = ComponentEdit()
        namespace_orig = {
            comp: {
                "filters": {
                    entity: OptionalFilter for entity in conf.get("filters", {})
                }
            }
            for comp, conf in config.items()
        }
        namespace = self.make_namespace(namespace_orig)  # type: ignore
        comp_edit.update_cli_namespace(namespace, {"pybids_inputs": config})
        for conf in config.values():
            assert conf.get("filters", {}) == {}

    @given(configs=two_configs(filters=False, paths=False))
    def test_wildcards_are_appended(self, configs: tuple[InputsConfig, InputsConfig]):
        comp_edit = ComponentEdit()
        config = configs[0]
        config_orig = copy.deepcopy(config)
        namespace_orig = configs[1]
        namespace = self.make_namespace(configs[1])
        comp_edit.update_cli_namespace(namespace, {"pybids_inputs": config})
        for comp, conf in config.items():
            assert conf.get("wildcards", []) == config_orig[comp].get(
                "wildcards", []
            ) + namespace_orig.get(comp, {}).get("wildcards", [])

    @given(configs=two_configs(filters=False, wildcards=False))
    def test_path_is_overwritten_and_resolved(
        self, configs: tuple[InputsConfig, InputsConfig]
    ):
        comp_edit = ComponentEdit()
        config = configs[0]
        config_orig = copy.deepcopy(config)
        namespace_orig = configs[1]
        namespace = self.make_namespace(configs[1])
        comp_edit.update_cli_namespace(namespace, {"pybids_inputs": config})
        for comp, conf in config.items():
            if (
                ref_path := namespace_orig.get(comp, {}).get("custom_path")
            ) is not None:
                custom_path = conf["custom_path"]  # type: ignore
                assert custom_path == Path(ref_path).resolve()
            else:
                assert conf.get("custom_path") == config_orig[comp].get("custom_path")

    @given(configs=two_configs())
    def test_namespace_is_cleared(self, configs: tuple[InputsConfig, InputsConfig]):
        comp_edit = ComponentEdit()
        config = configs[0]
        namespace = self.make_namespace(configs[1])
        comp_edit.update_cli_namespace(namespace, {"pybids_inputs": config})
        assert namespace == {}
