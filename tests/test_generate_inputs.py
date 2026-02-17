# ruff: noqa: PLR2004
from __future__ import annotations

import filecmp
import functools as ft
import itertools as it
import keyword
import logging
import os
import shutil
import sys
import tempfile
import warnings
from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path, PosixPath
from typing import Any, NamedTuple, TypeVar

import attrs
import more_itertools as itx
import pytest
from bids.layout import BIDSLayout
from hypothesis import HealthCheck, assume, example, given, settings
from hypothesis import strategies as st
from pytest_mock import MockerFixture

from snakebids.core import input_generation
from snakebids.core._querying import PostFilter, UnifiedFilter, get_matching_files
from snakebids.core.datasets import BidsComponent, BidsDataset
from snakebids.core.input_generation import (
    _all_custom_paths,
    _gen_bids_layout,
    _get_components,
    _is_local_relative,
    _normalize_database_args,
    _parse_bids_path,
    _parse_custom_path,
    generate_inputs,
)
from snakebids.exceptions import ConfigError, PybidsError, RunError
from snakebids.paths._presets import bids
from snakebids.snakemake_compat import expand as sb_expand
from snakebids.types import InputConfig, InputsConfig
from snakebids.utils.containers import MultiSelectDict
from snakebids.utils.utils import (
    DEPRECATION_FLAG,
    BidsEntity,
    BidsParseError,
    zip_list_eq,
)
from tests import strategies as sb_st
from tests.helpers import (
    Benchmark,
    BidsListCompare,
    allow_function_scoped,
    component_via_product,
    create_dataset,
    create_snakebids_config,
    get_bids_path,
    get_zip_list,
    reindex_dataset,
)

T = TypeVar("T")


def _not_deprecated(s: str):
    return not (s.startswith(DEPRECATION_FLAG) and s.endswith(DEPRECATION_FLAG))


class TestNormalizeDatabaseArgs:
    @given(
        pybidsdb_dir=st.text().filter(_not_deprecated) | st.none(),
        pybidsdb_reset=st.booleans(),
    )
    def test_normal_calls_give_no_warnings(
        self, pybidsdb_dir: str | None, pybidsdb_reset: bool
    ):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            _normalize_database_args(pybidsdb_dir, pybidsdb_reset, None, None)

    @given(
        pybidsdb_dir=st.text().filter(_not_deprecated) | st.none(),
        pybidsdb_reset=st.booleans(),
        pybids_database_dir=st.text().filter(_not_deprecated),
    )
    def test_old_dir_param_gives_warning(
        self, pybidsdb_dir: str | None, pybidsdb_reset: bool, pybids_database_dir: str
    ):
        with pytest.warns(UserWarning, match="`pybids_database_dir`"):
            _normalize_database_args(
                pybidsdb_dir, pybidsdb_reset, pybids_database_dir, None
            )

    @given(
        pybidsdb_dir=st.text().filter(_not_deprecated) | st.none(),
        pybidsdb_reset=st.booleans(),
        pybids_database_reset=st.booleans(),
    )
    def test_old_reset_param_gives_warning(
        self,
        pybidsdb_dir: str | None,
        pybidsdb_reset: bool,
        pybids_database_reset: bool,
    ):
        with pytest.warns(UserWarning, match="`pybids_reset_database`"):
            _normalize_database_args(
                pybidsdb_dir, pybidsdb_reset, None, pybids_database_reset
            )

    @given(
        pybids_database_dir=st.text().filter(_not_deprecated),
        pybids_database_reset=st.booleans(),
    )
    def test_old_params_passed_on_to_new(
        self,
        pybids_database_dir: str,
        pybids_database_reset: bool,
    ):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            assert (
                pybids_database_dir,
                pybids_database_reset,
            ) == _normalize_database_args(
                None, None, pybids_database_dir, pybids_database_reset
            )

    @given(
        pybidsdb_reset=st.booleans() | st.none(),
        pybids_database_reset=st.booleans() | st.none(),
    )
    def test_second_return_never_none(
        self,
        pybidsdb_reset: bool | None,
        pybids_database_reset: bool | None,
    ):
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=r"The parameter `pybids_reset_database` in generate_inputs\(\) "
                "is deprecated and will be removed in the next release. To reset the "
                "pybids database, use the `pybidsdb_reset` parameter instead",
                category=UserWarning,
            )
            assert (
                _normalize_database_args(
                    None,
                    pybidsdb_reset,
                    None,
                    pybids_database_reset,
                )[1]
                is not None
            )

    @given(pybidsdb_dir=st.text().filter(_not_deprecated))
    def test_deprecated_dir_raises_warning(self, pybidsdb_dir: str):
        pybidsdb_dir_ = DEPRECATION_FLAG + pybidsdb_dir + DEPRECATION_FLAG
        with pytest.warns(UserWarning, match="`pybids_db_dir`"):
            assert (
                _normalize_database_args(pybidsdb_dir_, None, None, None)[0]
                == pybidsdb_dir
            )

    @given(pybidsdb_reset=st.booleans())
    def test_deprecated_reset_raises_warning(self, pybidsdb_reset: bool):
        pybidsdb_reset_ = DEPRECATION_FLAG + str(int(pybidsdb_reset)) + DEPRECATION_FLAG
        with pytest.warns(UserWarning, match="`pybids_db_reset`"):
            assert (
                _normalize_database_args(None, pybidsdb_reset_, None, None)[1]
                == pybidsdb_reset
            )

    @given(pybidsdb_reset=st.text().filter(_not_deprecated))
    def test_non_deprecated_text_in_reset_raises_error(self, pybidsdb_reset: bool):
        with pytest.raises(TypeError):
            _normalize_database_args(None, pybidsdb_reset, None, None)


def test_regex_search_removed_from_filters():
    assert not len(UnifiedFilter.from_filter_dict({"regex_search": "foo"}).prefilters)


@given(
    filters=st.dictionaries(st.text(), st.text() | st.booleans() | st.lists(st.text()))
)
def test_get_matching_files_skips_get_when_empty_prefilter(filters: dict[str, Any]):
    assert (
        get_matching_files(
            ...,  # type: ignore
            UnifiedFilter.from_filter_dict({**filters, "foo": []}),
        )
        == []
    )


def test_attribute_errors_from_pybids_qualified_and_raised():
    with pytest.raises(PybidsError, match="Pybids has encountered a problem"):
        get_matching_files(..., UnifiedFilter.from_filter_dict({}))  # type: ignore


class TestFilterBools:
    def disambiguate_components(self, dataset: BidsDataset):
        assert len(dataset) == 2
        comp1, comp2 = dataset.values()
        return sorted([comp1, comp2], key=lambda comp: len(comp.entities.keys()))

    def get_extra_entity(self, dataset: BidsDataset) -> str:
        def set_difference(set_: set[str], comp: BidsComponent) -> set[str]:
            return set_.symmetric_difference(comp.entities.keys())

        return itx.one(
            ft.reduce(
                set_difference,
                dataset.values(),
                set[str](),
            )
        )

    @settings(
        deadline=800,
        suppress_health_check=[
            HealthCheck.function_scoped_fixture,
            HealthCheck.too_slow,
        ],
    )
    @given(dataset=sb_st.datasets())
    def test_ambiguous_paths_with_extra_entities_leads_to_error(
        self, tmpdir: Path, dataset: BidsDataset
    ):
        root = tempfile.mkdtemp(dir=tmpdir)
        create_dataset(root, dataset)
        shorter, _ = self.disambiguate_components(dataset)
        pybids_inputs: InputsConfig = {
            shorter.input_name: {
                "wildcards": [
                    BidsEntity.from_tag(wildcard).entity
                    for wildcard in shorter.input_wildcards
                ],
            }
        }

        with pytest.raises(ConfigError):
            generate_inputs(root, pybids_inputs)

    @settings(
        deadline=800,
        suppress_health_check=[
            HealthCheck.function_scoped_fixture,
            HealthCheck.too_slow,
        ],
    )
    @given(dataset=sb_st.datasets())
    def test_ambiguous_paths_with_missing_entity_leads_to_error(
        self, tmpdir: Path, dataset: BidsDataset
    ):
        root = tempfile.mkdtemp(dir=tmpdir)
        create_dataset(root, dataset)
        _, longer = self.disambiguate_components(dataset)
        pybids_inputs: InputsConfig = {
            longer.input_name: {
                "wildcards": [
                    BidsEntity.from_tag(wildcard).entity
                    for wildcard in longer.input_wildcards
                ],
            }
        }

        with pytest.raises(ConfigError):
            generate_inputs(root, pybids_inputs)

    @settings(
        deadline=800,
        suppress_health_check=[
            HealthCheck.function_scoped_fixture,
            HealthCheck.too_slow,
        ],
    )
    @given(dataset=sb_st.datasets(unique=True))
    def test_entity_excluded_when_filter_false(
        self, tmpdir: Path, dataset: BidsDataset
    ):
        root = tempfile.mkdtemp(dir=tmpdir)
        create_dataset(root, dataset)
        shorter, _ = self.disambiguate_components(dataset)
        extra_entity = self.get_extra_entity(dataset)
        pybids_inputs: InputsConfig = {
            shorter.input_name: {
                "wildcards": [
                    BidsEntity.from_tag(wildcard).entity
                    for wildcard in shorter.input_wildcards
                ],
                "filters": {BidsEntity.from_tag(extra_entity).entity: False},
            }
        }

        expected = BidsDataset(
            {
                shorter.input_name: attrs.evolve(
                    shorter, path=os.path.join(root, shorter.input_path)
                )
            }
        )

        data = generate_inputs(root, pybids_inputs)
        assert data == expected

    @example(
        dataset=BidsDataset(
            {
                "0": BidsComponent(
                    name="0",
                    path="ce-{ce}_space-{space}",
                    zip_lists={"ce": ["0"], "space": ["0"]},
                ),
                "1": BidsComponent(name="1", path="ce-{ce}", zip_lists={"ce": ["0"]}),
            }
        )
    )
    @example(
        dataset=BidsDataset(
            {
                "1": BidsComponent(
                    name="1",
                    path="sub-{subject}/{datatype}/sub-{subject}",
                    zip_lists={"subject": ["0"], "datatype": ["anat"]},
                ),
                "0": BidsComponent(
                    name="0",
                    path="sub-{subject}/sub-{subject}",
                    zip_lists={"subject": ["0"]},
                ),
            }
        )
    )
    @example(
        dataset=BidsDataset(
            {
                "1": BidsComponent(
                    name="1",
                    path="sub-{subject}/sub-{subject}_{suffix}.foo",
                    zip_lists={"subject": ["0"], "suffix": ["bar"]},
                ),
                "0": BidsComponent(
                    name="0",
                    path="{suffix}.foo",
                    zip_lists={"suffix": ["bar"]},
                ),
            }
        )
    )
    @settings(
        deadline=800,
        suppress_health_check=[
            HealthCheck.function_scoped_fixture,
            HealthCheck.too_slow,
        ],
    )
    @given(dataset=sb_st.datasets(unique=True))
    def test_entity_excluded_when_filter_true(self, tmpdir: Path, dataset: BidsDataset):
        root = tempfile.mkdtemp(dir=tmpdir)
        create_dataset(root, dataset)
        _, longer = self.disambiguate_components(dataset)
        extra_entity = self.get_extra_entity(dataset)
        pybids_inputs: InputsConfig = {
            longer.input_name: {
                "wildcards": [
                    BidsEntity.from_tag(wildcard).entity
                    for wildcard in longer.input_wildcards
                ],
                "filters": {BidsEntity.from_tag(extra_entity).entity: True},
            }
        }

        expected = BidsDataset(
            {
                longer.input_name: attrs.evolve(
                    longer, path=os.path.join(root, longer.input_path)
                )
            }
        )

        data = generate_inputs(root, pybids_inputs)
        assert data == expected

    @allow_function_scoped
    @given(
        template=sb_st.bids_components(
            name="template", restrict_patterns=True, unique=True, extra_entities=False
        ),
        entity=sb_st.bids_entity(path_safe=True),
        data=st.data(),
    )
    def test_text_filter_selects_paths_when_in_list_with_false(
        self,
        tmpdir: Path,
        template: BidsComponent,
        entity: BidsEntity,
        data: st.DataObject,
    ):
        """Also test here that the filters can pick out a desired value from a decoy"""

        def add_entity(component: BidsComponent, entity: str, value: str):
            return get_bids_path(component.wildcards, **{entity: value})

        assume(entity.wildcard not in template.wildcards)
        root = tempfile.mkdtemp(dir=tmpdir)

        target = data.draw(sb_st.bids_value(entity.match))
        decoy = data.draw(sb_st.bids_value(entity.match))
        assume(target != decoy)
        # This is a brute-force way of dealing with the fact that values for `run` are
        # converted into int by pybids before querying. So just prevent the decoy from
        # ever looking like the same number as target
        if target.isdecimal() and decoy.isdecimal():
            assume(int(target) != int(decoy))
        dataset = BidsDataset.from_iterable(
            [
                attrs.evolve(
                    template,
                    name="target",
                    path=os.path.join(
                        root, add_entity(template, entity.wildcard, target)
                    ),
                ),
                attrs.evolve(
                    template,
                    name="decoy",
                    path=os.path.join(
                        root, add_entity(template, entity.wildcard, decoy)
                    ),
                ),
            ]
        )
        create_dataset("", dataset)
        pybids_inputs: InputsConfig = {
            "target": {
                "wildcards": [
                    BidsEntity.from_tag(wildcard).entity
                    for wildcard in template.wildcards
                ],
                "filters": {entity.entity: [target, False]},
            }
        }

        result = generate_inputs(root, pybids_inputs)
        assert result == BidsDataset({"target": dataset["target"]})

    @allow_function_scoped
    @given(
        template=sb_st.bids_components(
            name="template", restrict_patterns=True, unique=True, extra_entities=False
        ),
        entity=sb_st.bids_entity(path_safe=True),
        data=st.data(),
    )
    def test_filter_false_in_list_selects_paths(
        self,
        tmpdir: Path,
        template: BidsComponent,
        entity: BidsEntity,
        data: st.DataObject,
    ):
        def add_entity(component: BidsComponent, entity: str, value: str):
            return get_bids_path(component.wildcards, **{entity: value})

        assume(entity.wildcard not in template.wildcards)
        target = data.draw(sb_st.bids_value(entity.match))
        decoy = data.draw(sb_st.bids_value(entity.match))
        assume(target != decoy)
        # See note in test_filter_works_when_false_in_list
        if target.isdecimal() and decoy.isdecimal():
            assume(int(target) != int(decoy))

        root = tempfile.mkdtemp(dir=tmpdir)
        dataset = BidsDataset.from_iterable(
            [
                attrs.evolve(
                    template,
                    path=os.path.join(root, template.path),
                ),
                attrs.evolve(
                    template,
                    name="decoy",
                    path=os.path.join(
                        root, add_entity(template, entity.wildcard, decoy)
                    ),
                ),
            ]
        )
        create_dataset("", dataset)
        pybids_inputs: InputsConfig = {
            "template": {
                "wildcards": [
                    BidsEntity.from_tag(wildcard).entity
                    for wildcard in template.wildcards
                ],
                "filters": {entity.entity: [target, False]},
            }
        }

        result = generate_inputs(root, pybids_inputs)
        assert result == BidsDataset({"template": dataset["template"]})


class TestFilterMethods:
    @pytest.mark.parametrize(
        ("entities", "query"),
        [
            (
                {
                    "subject": ["query", "qber7", "Query", "queryof"],
                    "tracer": ["1", "2"],
                },
                ["q.er[y7]", "1"],
            ),
            (
                {
                    "subject": ["0", "00", "01"],
                    "mt": ["on"],
                },
                ["0+", "on"],
            ),
        ],
    )
    def test_regex_match_selects_paths(
        self, tmpdir: Path, entities: dict[str, list[str]], query: list[str]
    ):
        component = component_via_product(**entities)
        dataset = BidsDataset.from_iterable([component])
        create_dataset(tmpdir, dataset)
        pybids_inputs: InputsConfig = {
            "template": {
                "wildcards": list(entities),
                "filters": {key: {"match": q} for key, q in zip(entities, query)},  # noqa: B905
            }
        }
        result = generate_inputs(tmpdir, pybids_inputs)
        assert len(result["template"].expand()) == 2

    @pytest.mark.parametrize(
        "entities",
        [
            {"subject": ["query", "pquerya", "Query"], "tracer": ["1", "2"]},
            {"proc": ["on", "none"], "subject": ["0", "1"]},
        ],
    )
    def test_regex_search_selects_paths(
        self, tmpdir: Path, entities: dict[str, list[str]]
    ):
        component = component_via_product(**entities)
        dataset = BidsDataset.from_iterable([component])
        create_dataset(tmpdir, dataset)
        pybids_inputs: InputsConfig = {
            "template": {
                "wildcards": list(entities),
                "filters": {key: {"search": val[0]} for key, val in entities.items()},
            }
        }
        result = generate_inputs(tmpdir, pybids_inputs)
        assert len(result["template"].expand()) == 2

    @pytest.mark.parametrize(
        "entities",
        [
            {"subject": ["001", "002"], "acquisition": ["1", "11", "12", "21"]},
            {"subject": ["0", "00"], "mt": ["on", "off"]},
        ],
    )
    def test_get_method_selects_via_direct_matching(
        self, tmpdir: Path, entities: dict[str, list[str]]
    ):
        component = component_via_product(**entities)
        dataset = BidsDataset.from_iterable([component])
        create_dataset(tmpdir, dataset)
        pybids_inputs: InputsConfig = {
            "template": {
                "filters": {key: {"get": val[0]} for key, val in entities.items()},
            }
        }
        result = generate_inputs(tmpdir, pybids_inputs)
        assert len(result["template"].expand()) == 1

    def test_combining_match_and_get_selects_correct_paths(self, tmpdir: Path):
        component = component_via_product(
            subject=["001", "002"], acq=["1", "11", "12", "21"]
        )
        dataset = BidsDataset.from_iterable([component])
        create_dataset(tmpdir, dataset)
        pybids_inputs: InputsConfig = {
            "template": {
                "filters": {
                    "subject": {"get": "001"},
                    "acquisition": {"match": "\\d"},
                },
            }
        }
        result = generate_inputs(tmpdir, pybids_inputs)
        assert len(result["template"].expand()) == 1

    @given(
        methods=st.lists(
            st.sampled_from(["match", "get", "search"]), unique=True, min_size=2
        )
    )
    def test_filter_with_multiple_methods_raises_error(self, methods: list[str]):
        pybids_inputs: InputsConfig = {
            "template": {
                "filters": {
                    "foo": dict.fromkeys(methods, "foo")  # type: ignore
                },
            }
        }
        with pytest.raises(ConfigError, match="may not have more than one key"):
            generate_inputs("scripts", pybids_inputs)

    def test_filter_with_no_methods_raises_error(self, tmpdir: Path):
        pybids_inputs: InputsConfig = {
            "template": {
                "filters": {"foo": {}},
            }
        }
        with pytest.raises(ConfigError, match="was not given any keys"):
            generate_inputs(tmpdir, pybids_inputs)

    @given(method=st.text().filter(lambda s: s not in {"get", "match", "search"}))
    def test_filter_with_invalid_method_raises_error(self, method: str):
        pybids_inputs: InputsConfig = {
            "template": {
                "filters": {"foo": {method: []}},  # type: ignore
            }
        }
        with pytest.raises(ConfigError, match="Invalid query method specified"):
            generate_inputs("scripts", pybids_inputs)


@settings(
    deadline=800,
    suppress_health_check=[
        HealthCheck.function_scoped_fixture,
        HealthCheck.too_slow,
    ],
    max_examples=1,
)
@given(dataset=sb_st.datasets_one_comp(unique=True))
def test_duplicate_wildcards_does_not_create_error(dataset: BidsDataset, tmpdir: Path):
    root = tempfile.mkdtemp(dir=tmpdir)
    rooted = BidsDataset.from_iterable(
        attrs.evolve(comp, path=os.path.join(root, comp.path))
        for comp in dataset.values()
    )
    create_dataset(Path("/"), rooted)
    config = create_snakebids_config(dataset)
    wildcards = itx.first(config.values()).get("wildcards", [])
    wildcards.append(wildcards[0])
    reindexed = generate_inputs(
        root,
        config,
    )
    assert reindexed == rooted
    assert reindexed.layout is not None


class TestAbsentConfigEntries:
    def get_entities(self, root: Path):
        # Generate directory
        entities = {"subject": ["001", "002"], "acq": sorted(["foo", "bar"])}
        zip_list: dict[str, list[str]] = defaultdict(list)
        for e in it.product(*entities.values()):
            d = dict(zip(entities.keys(), e, strict=True))
            for key, val in d.items():
                zip_list[key].append(val)
            path = Path(bids(root, datatype="anat", suffix="T1w.nii.gz", **d))
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()
        return entities, MultiSelectDict(zip_list)

    def test_missing_filters(self, tmpdir: Path):
        _, zip_list = self.get_entities(tmpdir)

        # create config
        derivatives = False
        pybids_inputs: InputsConfig = {
            "t1": {
                "wildcards": ["acquisition", "subject"],
            }
        }

        # Simplest case -- one input type, using pybids
        config = generate_inputs(
            pybids_inputs=pybids_inputs,
            bids_dir=tmpdir,
            derivatives=derivatives,
            pybids_config=str(Path(__file__).parent / "data" / "custom_config.json"),
        )
        template = BidsDataset(
            {"t1": BidsComponent(name="t1", path=config["t1"].path, zip_lists=zip_list)}
        )
        # Order of the subjects is not deterministic
        assert template == config
        assert config.subj_wildcards == {"subject": "{subject}"}

    def test_missing_wildcards(self, tmpdir: Path):
        self.get_entities(tmpdir)

        # create config
        derivatives = False
        pybids_inputs: InputsConfig = {
            "t1": {
                "filters": {"acquisition": "foo", "subject": "001"},
            }
        }

        # Simplest case -- one input type, using pybids
        config = generate_inputs(
            bids_dir=tmpdir,
            pybids_inputs=pybids_inputs,
            derivatives=derivatives,
            pybids_config=str(Path(__file__).parent / "data" / "custom_config.json"),
        )
        template = BidsDataset(
            {"t1": BidsComponent(name="t1", path=config["t1"].path, zip_lists={})}
        )
        assert template == config
        assert config.subj_wildcards == {"subject": "{subject}"}


class PathEntities(NamedTuple):
    entities: dict[str, list[str]]
    template: Path
    filters: dict[str, list[str]]


@st.composite
def path_entities(draw: st.DrawFn):
    """Generate path wildcard entities and corresponding values

    This has three main outputs. The first is a dict of wildcards->values, where values
    are sample values of wildcards. Wildcards are substituteable, {brace-enclosed} str
    for use in string.format.

    The wildcards in the entity dict produced above are combined to form a path name:

    wildcard1-{wildcard1}_wildcard2-{wildcard2}_...

    In addition to this name, a random subsample of wildcards will be selected as dir
    entities. These, along with the path name, will be combined to return a path
    template:

    {wildcard2}/{wildcard5}/wildcard1-{wildcard1}_wildcard2-{wildcard2}_...

    Finally, a subsample of wildcards, along with a subset of their respective values,
    will be selected as filters and returned as a dict of wildcard->values.

    Parameters
    ----------
    draw : st.DrawFn
        [description]

    Returns
    -------
    [type]
        [description]
    """
    # TODO: Remove restrictions on characters to accurately reflect what the function
    #       should allow
    valid_chars = st.characters(
        min_codepoint=48, max_codepoint=122, whitelist_categories=["Ll", "Lu"]
    )
    # We need to explicitly exclude keywords here because the current implementation
    # of glob_wildcards uses a named tuple, which doesn't allow keywords as attributes.
    path_text = st.text(valid_chars, min_size=1).filter(
        lambda s: not keyword.iskeyword(s)
    )
    entities = draw(
        st.dictionaries(
            path_text,
            st.lists(path_text, min_size=1, max_size=2, unique=True),
            min_size=1,
            max_size=5,
        )
    )

    def get_subset(of: Iterable[T]) -> list[T]:
        return draw(
            st.lists(st.sampled_from([*of]), unique=True, max_size=len(entities))
        )

    # Collect dir_entities and filters
    dir_entities = get_subset(entities.keys())
    filtered_entities = get_subset(entities.keys())
    filter_selections = [get_subset(entities[entity]) for entity in filtered_entities]
    filters = dict(zip(filtered_entities, filter_selections, strict=True))

    # Compose the path template
    dir_template = Path(*(f"{{{entity}}}" for entity in dir_entities))
    name_template = "_".join(f"{entity}-{{{entity}}}" for entity in entities)
    template = dir_template / name_template

    return PathEntities(entities, template, filters)


class TestCustomPaths:
    @pytest.fixture
    def temp_dir(self, fakefs_tmpdir: Path, bids_fs: Path):
        return fakefs_tmpdir

    def generate_test_directory(
        self, entities: dict[str, list[str]], template: Path, tmpdir: Path
    ):
        root = Path(tempfile.mkdtemp(prefix="hypothesis-", dir=tmpdir))
        # Generate fake directory structure
        for values in it.product(*entities.values()):
            name_value = dict(zip(entities.keys(), values, strict=True))
            path = str(template).format(**name_value)
            (root / path).parent.mkdir(parents=True, exist_ok=True)
            (root / path).touch()
        return root / template

    T = TypeVar("T")

    def test_benchmark_test_custom_paths(self, benchmark: Benchmark, tmp_path: Path):
        entities = {"A": ["A", "B", "C"], "B": ["1", "2", "3"]}
        template = Path("{A}/A-{A}_B-{B}")
        test_path = self.generate_test_directory(entities, template, tmp_path)
        benchmark(_parse_custom_path, test_path, UnifiedFilter.from_filter_dict({}))

    @allow_function_scoped
    @given(path_entities=path_entities())
    def test_collects_all_paths_when_no_filters(
        self,
        path_entities: PathEntities,
        temp_dir: Path,
    ):
        entities, template, _ = path_entities
        test_path = self.generate_test_directory(entities, template, temp_dir)

        # Test without any filters
        result = _parse_custom_path(test_path, UnifiedFilter.from_filter_dict({}))
        zip_lists = get_zip_list(entities, it.product(*entities.values()))
        assert BidsComponent(
            name="foo", path=get_bids_path(zip_lists), zip_lists=zip_lists
        ) == BidsComponent(
            name="foo", path=get_bids_path(result), zip_lists=MultiSelectDict(result)
        )

    @settings(deadline=400, suppress_health_check=[HealthCheck.function_scoped_fixture])
    @given(path_entities=path_entities())
    def test_collects_only_filtered_entities(
        self,
        path_entities: PathEntities,
        temp_dir: Path,
    ):
        entities, template, filters = path_entities
        test_path = self.generate_test_directory(entities, template, temp_dir)

        # Test with filters
        result_filtered = MultiSelectDict(
            _parse_custom_path(test_path, UnifiedFilter.from_filter_dict(filters))
        )
        zip_lists = MultiSelectDict(
            # Start with empty lists for each key, otherwise keys will be missing
            {key: list[str]() for key in entities}
            |
            # Override entities with relevant filters before making zip lists
            get_zip_list(entities, it.product(*(entities | filters).values()))
        )
        assert BidsComponent(
            name="foo", path=get_bids_path(zip_lists), zip_lists=zip_lists
        ) == BidsComponent(
            name="foo", path=get_bids_path(result_filtered), zip_lists=result_filtered
        )

    @settings(
        deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
    )
    @given(path_entities=path_entities())
    def test_collect_all_but_filters_when_exclusion_filters_used(
        self,
        path_entities: PathEntities,
        temp_dir: Path,
    ):
        entities, template, filters = path_entities
        test_path = self.generate_test_directory(entities, template, temp_dir)
        # Test with exclusion filters
        exclude_filters = PostFilter()
        for key, values in filters.items():
            exclude_filters.add_filter(key, None, values)
        result_excluded = _parse_custom_path(
            test_path, UnifiedFilter.from_filter_dict({}, exclude_filters)
        )

        entities_excluded = {
            entity: [value for value in values if value not in filters.get(entity, [])]
            for entity, values in entities.items()
        }
        zip_lists = MultiSelectDict(
            # Start with empty lists for each key, otherwise keys will be missing
            {key: list[str]() for key in entities}
            | get_zip_list(entities, it.product(*entities_excluded.values())),
        )

        assert zip_list_eq(zip_lists, result_excluded)

    @given(
        boolean=st.booleans(),
        filter=st.none() | st.lists(st.text()),
        path_entities=path_entities(),
    )
    @allow_function_scoped
    def test_errors_when_bools_given_as_filters(
        self,
        temp_dir: Path,
        path_entities: PathEntities,
        boolean: bool,
        filter: list[str] | None,
    ):
        entities, template, _ = path_entities
        test_path = self.generate_test_directory(entities, template, temp_dir)
        with pytest.raises(ValueError, match="Boolean filters in items with custom "):
            _parse_custom_path(
                test_path,
                UnifiedFilter.from_filter_dict(
                    {"foo": boolean if filter is None else [*filter, boolean]}
                ),
            )


def test_custom_pybids_config(tmpdir: Path):
    # Generate directory
    for i in range(2):
        path = Path(
            bids(
                tmpdir, datatype="anat", subject="001", foo=str(i), suffix="T1w.nii.gz"
            )
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()

    # create config
    derivatives = False
    pybids_inputs: InputsConfig = {
        "t1": {
            "filters": {"suffix": "T1w"},
            "wildcards": ["acquisition", "subject", "foo"],
        }
    }

    # Simplest case -- one input type, using pybids
    result = generate_inputs(
        pybids_inputs=pybids_inputs,
        bids_dir=tmpdir,
        derivatives=derivatives,
        pybids_config=str(Path(__file__).parent / "data" / "custom_config.json"),
    )
    template = BidsDataset(
        {
            "t1": BidsComponent(
                name="t1",
                path=bids(
                    tmpdir,
                    datatype="anat",
                    subject="{subject}",
                    foo="{foo}",
                    suffix="T1w.nii.gz",
                ),
                zip_lists={"foo": ["0", "1"], "subject": ["001", "001"]},
            )
        }
    )
    assert template == result
    assert result["t1"].wildcards == {"foo": "{foo}", "subject": "{subject}"}
    # Order of the subjects is not deterministic
    assert result.subj_wildcards == {"subject": "{subject}"}


def test_nonstandard_custom_pybids_config(tmpdir: Path):
    # Generate directory
    for i in range(2):
        path = Path(
            bids(
                tmpdir, datatype="anat", subject="001", foo=str(i), suffix="T1w.nii.gz"
            )
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()

    # create config
    derivatives = False
    pybids_inputs: InputsConfig = {
        "t1": {
            "filters": {"suffix": "T1w"},
            "wildcards": ["acquisition", "subject", "foobar"],
        }
    }

    # Simplest case -- one input type, using pybids
    with pytest.raises(ConfigError):
        generate_inputs(
            pybids_inputs=pybids_inputs,
            bids_dir=tmpdir,
            derivatives=derivatives,
            pybids_config=(
                str(Path(__file__).parent / "data" / "custom_config_nonstandard.json")
            ),
        )


def test_index_metadata(mocker: MockerFixture):
    spy = mocker.spy(input_generation, "BIDSLayoutIndexer")
    mocker.patch.object(input_generation, "BIDSLayout", side_effect=ValueError)

    # Simplest case -- one input type, using pybids
    with pytest.raises(ValueError):  # noqa
        generate_inputs(
            pybids_inputs={"foo": {}},
            bids_dir=...,  # type: ignore
            derivatives=...,  # type: ignore
            index_metadata=True,
        )

    spy.assert_called_once_with(
        validate=False,
        index_metadata=True,
    )


def test_t1w():
    # create config
    real_bids_dir = "tests/data/bids_t1w"
    derivatives = False
    pybids_inputs: InputsConfig = {
        "t1": {
            "filters": {"suffix": "T1w"},
            "wildcards": ["acquisition", "subject", "session", "run"],
        }
    }

    # Simplest case -- one input type, using pybids
    result = generate_inputs(
        pybids_inputs=pybids_inputs,
        bids_dir=real_bids_dir,
        derivatives=derivatives,
    )
    template = BidsDataset(
        {
            "t1": BidsComponent(
                name="t1",
                path=result["t1"].path,
                zip_lists={"acq": ["mprage", "mprage"], "subject": ["001", "002"]},
            )
        }
    )
    assert template == result

    # Order of the subjects is not deterministic
    assert result.subjects in [["001", "002"], ["002", "001"]]
    assert result.sessions == []
    assert result.subj_wildcards == {"subject": "{subject}"}

    pybids_inputs_suffix: InputsConfig = {
        "scan": {
            "filters": {},
            "wildcards": [
                "acquisition",
                "subject",
                "session",
                "run",
                "suffix",
            ],
        }
    }
    result = generate_inputs(
        pybids_inputs=pybids_inputs_suffix,
        bids_dir=real_bids_dir,
        derivatives=derivatives,
        participant_label="001",
    )
    assert result["scan"].entities == {
        "acq": ["mprage"],
        "subject": ["001"],
        "suffix": ["T1w"],
    }
    template = BidsDataset(
        {
            "scan": BidsComponent(
                name="scan",
                path=result["scan"].path,
                zip_lists={
                    "acq": [
                        "mprage",
                    ],
                    "subject": [
                        "001",
                    ],
                    "suffix": [
                        "T1w",
                    ],
                },
            )
        }
    )
    assert template == result

    assert result.subjects == ["001"]
    assert result.sessions == []
    assert result.subj_wildcards == {"subject": "{subject}"}

    # Two input types, specified by pybids or path override
    wildcard_path_t1 = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "data/bids_t1w",
        "sub-{subject}/anat/sub-{subject}_acq-{acq}_T1w.nii.gz",
    )
    wildcard_path_t2 = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "data/bids_t1w",
        "sub-{subject}/anat/sub-{subject}_T2w.nii.gz",
    )
    pybids_inputs = {
        "t1": {
            "filters": {"suffix": "T1w"},
            "wildcards": ["acquisition", "subject", "session", "run"],
        },
        "t2": {
            "filters": {"suffix": "T2w"},
            "wildcards": ["acquisition", "subject", "session", "run"],
        },
    }
    bids_dir = real_bids_dir

    # Want to test both inputs from layout, both inputs from custom path, and
    # one of each. This setup should produce the same results every time.
    for idx in range(4):
        if idx == 1:
            pybids_inputs["t1"]["custom_path"] = wildcard_path_t1
        elif idx == 2:
            pybids_inputs["t2"]["custom_path"] = wildcard_path_t2
        elif idx == 3:
            pybids_inputs["t1"]["custom_path"] = wildcard_path_t1
            pybids_inputs["t2"]["custom_path"] = wildcard_path_t2
            # TODO: Allow arbitrary paths to work when all custom paths are used
            # bids_dir = "-"
        result = generate_inputs(
            pybids_inputs=pybids_inputs,
            bids_dir=bids_dir,
            derivatives=derivatives,
        )
        template = BidsDataset(
            {
                "t1": BidsComponent(
                    name="t1",
                    path=result["t1"].path,
                    zip_lists={
                        "acq": ["mprage", "mprage"],
                        "subject": ["001", "002"],
                    },
                ),
                "t2": BidsComponent(
                    name="t2",
                    path=result["t2"].path,
                    zip_lists={"subject": ["002"]},
                ),
            }
        )
        # Order of the subjects is not deterministic
        assert result.subjects in [["001", "002"], ["002", "001"]]

        assert result.sessions == []
        assert result.subj_wildcards == {"subject": "{subject}"}


def test_t1w_with_dict():
    # create config
    real_bids_dir = "tests/data/bids_t1w"
    derivatives = False
    pybids_inputs: InputsConfig = {
        "t1": {
            "filters": {"suffix": "T1w"},
            "wildcards": ["acquisition", "subject", "session", "run"],
        }
    }

    # Simplest case -- one input type, using pybids
    config = generate_inputs(
        pybids_inputs=pybids_inputs,
        bids_dir=real_bids_dir,
        derivatives=derivatives,
        use_bids_inputs=False,
    )
    # Order of the subjects is not deterministic
    assert config["input_lists"] == BidsListCompare(
        {"t1": {"acq": ["mprage"], "subject": ["002", "001"]}}
    )
    assert config["input_zip_lists"] == {
        "t1": {"acq": ["mprage", "mprage"], "subject": ["001", "002"]}
    }
    assert config["input_wildcards"] == {"t1": {"acq": "{acq}", "subject": "{subject}"}}
    # Order of the subjects is not deterministic
    assert set(config["subjects"]) == {"002", "001"}
    assert config["sessions"] == []
    assert config["subj_wildcards"] == {"subject": "{subject}"}

    pybids_inputs_suffix: InputsConfig = {
        "scan": {
            "filters": {},
            "wildcards": [
                "acquisition",
                "subject",
                "session",
                "run",
                "suffix",
            ],
        }
    }
    config = generate_inputs(
        pybids_inputs=pybids_inputs_suffix,
        bids_dir=real_bids_dir,
        derivatives=derivatives,
        participant_label="001",
        use_bids_inputs=False,
    )
    assert config["input_lists"] == {
        "scan": {"acq": ["mprage"], "subject": ["001"], "suffix": ["T1w"]}
    }
    assert config["input_zip_lists"] == {
        "scan": {"acq": ["mprage"], "subject": ["001"], "suffix": ["T1w"]}
    }
    assert config["input_wildcards"] == {
        "scan": {"acq": "{acq}", "subject": "{subject}", "suffix": "{suffix}"}
    }
    assert config["subjects"] == ["001"]
    assert config["sessions"] == []
    assert config["subj_wildcards"] == {"subject": "{subject}"}


def test_get_lists_from_bids_raises_run_error():
    bids_layout = None
    pybids_inputs: InputsConfig = {
        "t1": {
            "filters": {"suffix": "T1w"},
            "wildcards": ["acquisition", "subject", "session", "run"],
        }
    }
    with pytest.raises(RunError):
        next(
            _get_components(
                bids_layout=bids_layout,
                inputs_config=pybids_inputs,
                postfilters=PostFilter(),
            )
        )


def test_get_lists_from_bids():
    bids_dir = "tests/data/bids_t1w"
    wildcard_path_t1 = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "data/bids_t1w",
        "sub-{subject}/anat/sub-{subject}_acq-{acq}_T1w.nii.gz",
    )
    wildcard_path_t2 = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "data/bids_t1w",
        "sub-{subject}/anat/sub-{subject}_T2w.nii.gz",
    )
    print(wildcard_path_t1)
    layout = BIDSLayout(bids_dir, validate=False)
    pybids_inputs: InputsConfig = {
        "t1": {
            "filters": {"suffix": "T1w"},
            "wildcards": ["acquisition", "subject", "session", "run"],
        },
        "t2": {
            "filters": {"suffix": "T2w"},
            "wildcards": ["acquisition", "subject", "session", "run"],
        },
    }

    # Want to test both inputs from layout, both inputs from custom path, and
    # one of each. This setup should produce the same results every time.
    for idx in range(4):
        if idx == 1:
            pybids_inputs["t1"]["custom_path"] = wildcard_path_t1
        elif idx == 2:
            pybids_inputs["t2"]["custom_path"] = wildcard_path_t2
        elif idx == 3:
            pybids_inputs["t1"]["custom_path"] = wildcard_path_t1
            pybids_inputs["t2"]["custom_path"] = wildcard_path_t2

        result = _get_components(
            bids_layout=layout, inputs_config=pybids_inputs, postfilters=PostFilter()
        )
        for bids_lists in result:
            if bids_lists.input_name == "t1":
                template = BidsComponent(
                    name="t1",
                    path=wildcard_path_t1,
                    zip_lists={
                        "acq": ["mprage", "mprage"],
                        "subject": ["001", "002"],
                    },
                )
                assert template == bids_lists
            elif bids_lists.input_name == "t2":
                assert bids_lists.input_path == wildcard_path_t2
                template = BidsComponent(
                    name="t2",
                    path=wildcard_path_t2,
                    zip_lists={
                        "subject": ["002"],
                    },
                )
                assert template == bids_lists


class TestGenBidsLayout:
    @settings(
        max_examples=1, suppress_health_check=[HealthCheck.function_scoped_fixture]
    )
    @given(dataset=sb_st.datasets())
    def test_gen_layout_returns_valid_dataset(self, dataset: BidsDataset, tmpdir: Path):
        create_dataset(tmpdir, dataset)
        assert _gen_bids_layout(
            bids_dir=tmpdir,
            derivatives=False,
            pybidsdb_dir=None,
            pybidsdb_reset=False,
            pybids_config=None,
        )

    def test_invalid_path_raises_error(self, tmpdir: Path):
        with pytest.raises(ValueError, match="BIDS root does not exist"):
            _gen_bids_layout(
                bids_dir=tmpdir / "foo",
                derivatives=False,
                pybidsdb_dir=None,
                pybidsdb_reset=False,
            )


@pytest.mark.parametrize("count", tuple(range(6)))
def test_all_custom_paths(count: int):
    config: InputsConfig = {}
    for i in range(count):
        config[str(i)] = {"wildcards": [], "filters": {}, "custom_path": "foo"}
    for i in range(count, 5):
        config[str(i)] = {
            "wildcards": [],
            "filters": {},
        }
    if count == 5:
        assert _all_custom_paths(config)
    else:
        assert not _all_custom_paths(config)


class TestRecogPathSchemes:
    PATH_AND_TYPES = (
        ("file", "RELATIVE"),
        ("hello", "RELATIVE"),
        ("gs", "RELATIVE"),
        ("./hello/world", "RELATIVE"),
        ("hello/world", "RELATIVE"),
        ("/hello/world", "ABSOLUTE"),
        ("gs://some/google/cloud/bucket", "NETWORK"),
        ("s3://some/s3/bucket", "NETWORK"),
    )

    @pytest.mark.parametrize(("path", "path_type"), PATH_AND_TYPES)
    def test_is_local_relative(self, path: str, path_type: str):
        isnet = path_type == "NETWORK"
        is_local_relative = path_type == "RELATIVE"

        # test the path itself, and the corresponding Path(path)
        assert _is_local_relative(path) == is_local_relative, (
            f"Path {path} fails is local relative path test."
        )
        if not isnet:
            assert _is_local_relative(Path(path)) == is_local_relative

    @pytest.mark.skipif(
        sys.version_info < (3, 12), reason="Path class has no with_segments()"
    )
    @pytest.mark.parametrize(
        ("path", "path_type"), [tup for tup in PATH_AND_TYPES if tup[1] == "RELATIVE"]
    )
    def test_path_subclassing(self, path: str, path_type: str):
        # Google cloud is not posix, for mocking purpose however we just
        # need a class that is a subclass of Path
        class MockGCSPath(PosixPath):
            def __init__(self, *pathsegments: str):  # pyright: ignore[reportInconsistentConstructor]
                super().__init__(*pathsegments)

            def __str__(self):  # __fspath__ calls __str__ by default
                return f"gs://{super().__str__()}"

        assert not _is_local_relative(MockGCSPath(path))


@example(
    dataset=BidsDataset(
        {
            "1": BidsComponent(
                name="1",
                path="sub-{subject}/sub-{subject}_{suffix}{extension}",
                zip_lists={
                    "subject": ["0"],
                    "suffix": ["0"],
                    "extension": [".0"],
                },
            ),
            "0": BidsComponent(
                name="0",
                path="sub-{subject}/sub-{subject}{extension}",
                zip_lists={
                    "subject": ["0"],
                    "extension": [".0"],
                },
            ),
        }
    ),
)
@settings(
    deadline=800,
    suppress_health_check=[
        HealthCheck.function_scoped_fixture,
        HealthCheck.too_slow,
    ],
)
@given(dataset=sb_st.datasets(unique=True))
def test_generate_inputs(dataset: BidsDataset, tmpdir: Path):
    root = tempfile.mkdtemp(dir=tmpdir)
    rooted = BidsDataset.from_iterable(
        attrs.evolve(comp, path=os.path.join(root, comp.path))
        for comp in dataset.values()
    )
    reindexed = reindex_dataset(root, rooted)
    assert reindexed == rooted
    assert reindexed.layout is not None


class _FakeBIDSLayout:
    def __init__(self, *args: Any, **kwargs: Any):
        pass

    def get(self, regex_search: bool, **kwargs: Any) -> list[str]:
        return []

    def __eq__(self, other: object):
        return isinstance(other, self.__class__)


class TestParticipantFiltering:
    @pytest.mark.parametrize(
        ("include", "exclude"),
        [
            ("01", None),
            (None, "01"),
            ("01", "02"),
            (["01", "02"], "02"),
            (None, ["02", "03"]),
            (["01", "02"], ["02", "03"]),
            (["01", "02"], []),
            ([], ["02", "03"]),
            ([], []),
        ],
    )
    def test_exclude_and_participant_label_filter_correctly(
        self,
        include: str | list[str] | None,
        exclude: str | list[str] | None,
        tmpdir: Path,
    ):
        zip_lists = {
            "subject": ["01", "02", "03"],
            "acq": ["x", "y", "z"],
        }
        component = BidsComponent(
            name="0", path=str(tmpdir / get_bids_path(zip_lists)), zip_lists=zip_lists
        )
        dataset = BidsDataset.from_iterable([component])

        reindexed = reindex_dataset(
            str(tmpdir),
            dataset,
            participant_label=include,
            exclude_participant_label=exclude,
        )
        expected = set(component.entities["subject"])
        expected -= set(itx.always_iterable(exclude))
        if include is not None:
            expected &= set(itx.always_iterable(include))
        assert set(itx.first(reindexed.values()).entities["subject"]) == expected

    @pytest.mark.parametrize("include", ["include", ["include"], None])
    @pytest.mark.parametrize("exclude", ["exclude", ["exclude"], None])
    @pytest.mark.parametrize("filters", ["filters", None])
    def test_participant_flags_and_filters_merged(
        self,
        include: str | list[str] | None,
        exclude: str | list[str] | None,
        filters: str | None,
        mocker: MockerFixture,
        caplog: pytest.LogCaptureFixture,
    ):
        mocker.patch.object(input_generation, "BIDSLayout", _FakeBIDSLayout)
        patch = mocker.patch.object(input_generation, "get_matching_files")
        component: InputConfig = {"filters": {"subject": filters}} if filters else {}
        with caplog.at_level(logging.ERROR, "snakebids.core.input_generation"):
            generate_inputs(
                "",
                {"": component},
                participant_label=include,
                exclude_participant_label=exclude,
            )
        pf = PostFilter()
        pf.add_filter("subject", include, exclude)
        uf = UnifiedFilter(component, pf)
        patch.assert_called_once_with(_FakeBIDSLayout(), uf)

    def test_exclude_participant_does_not_make_all_other_filters_regex(
        self, tmpdir: Path
    ):
        # Construct a small deterministic BIDS dataset with a subject entity and
        # at least one other entity so we can mutate that other entity's values.
        # Single-component dataset with a valid BIDS-style path and entities.
        mut_entity = "acq"
        zip_lists = {
            "subject": ["01", "02"],
            mut_entity: ["x", "y"],
        }
        component = BidsComponent(
            name="0", path=str(tmpdir / get_bids_path(zip_lists)), zip_lists=zip_lists
        )
        dataset = BidsDataset.from_iterable([component])

        # Create an extra set of paths by modifying one of the existing components to
        # put 'foo' after a set of entity values. If that filter gets changed to a
        # regex, all of the suffixed decoys will get picked up by pybids.
        zip_lists[mut_entity] = ["foo" + v for v in zip_lists[mut_entity]]
        for path in sb_expand(itx.first(dataset.values()).path, zip, **zip_lists):
            p = Path(path)
            p.parent.mkdir(parents=True, exist_ok=True)
            p.touch()

        reindexed = reindex_dataset(
            str(tmpdir), dataset, exclude_participant_label="01"
        )
        assert itx.first(reindexed.values()) == component.filter(subject="02")


# The content of the dataset is irrelevant to this test, so one example suffices
# but can't use extension because the custom path won't glob properly
@settings(max_examples=1, suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(dataset=sb_st.datasets_one_comp(blacklist_entities=["extension"], unique=True))
def test_when_all_custom_paths_no_layout_indexed(
    dataset: BidsDataset, tmpdir: Path, mocker: MockerFixture
):
    # Need to reset mocker at beginning because hypothesis may call this function
    # multiple times
    mocker.stopall()
    root = tempfile.mkdtemp(dir=tmpdir)
    rooted = BidsDataset.from_iterable(
        attrs.evolve(comp, path=os.path.join(root, comp.path))
        for comp in dataset.values()
    )

    spy = mocker.spy(BIDSLayout, "__init__")
    reindexed = reindex_dataset(root, rooted, use_custom_paths=True)
    assert reindexed == rooted
    assert reindexed.layout is None
    spy.assert_not_called()


class TestParseBidsPath:
    @given(
        component=sb_st.bids_components(max_values=1, restrict_patterns=True),
        scheme=sb_st.schemes() | st.none(),
    )
    def test_splits_wildcards_from_path(
        self, component: BidsComponent, scheme: str | None
    ):
        path = component.expand()[0]
        if scheme is not None:
            path = f"{scheme}{path}"
        entities = [BidsEntity.normalize(e).entity for e in component.zip_lists]
        tpl_path, matches = _parse_bids_path(path, entities)
        assert tpl_path.format(**matches) == path

    @given(
        component=sb_st.bids_components(max_values=1, restrict_patterns=True),
        scheme=sb_st.schemes() | st.none(),
    )
    def test_one_match_found_for_each_entity(
        self, component: BidsComponent, scheme: str | None
    ):
        path = component.expand()[0]
        if scheme is not None:
            path = f"{scheme}{path}"
        entities = [BidsEntity.normalize(e).entity for e in component.zip_lists]
        _, matches = _parse_bids_path(path, entities)
        assert set(matches.items()) == {
            (key, val[0]) for key, val in component.zip_lists.items()
        }

    @given(
        component=sb_st.bids_components(
            max_values=1, restrict_patterns=True, extra_entities=False
        ),
        scheme=sb_st.schemes() | st.none(),
        entity=sb_st.bids_entity(),
    )
    def test_missing_match_leads_to_error(
        self, component: BidsComponent, scheme: str | None, entity: BidsEntity
    ):
        path = component.expand()[0]
        if scheme is not None:
            path = f"{scheme}{path}"
        entities = [BidsEntity.normalize(e).entity for e in component.zip_lists]
        assume(entity.entity not in entities)
        with pytest.raises(BidsParseError) as err:
            _parse_bids_path(path, it.chain(entities, [entity.entity]))
        assert err.value.entity == entity


class TestDB:
    @pytest.fixture(autouse=True)
    def _start(self, tmp_path: Path):
        self.tmpdir = str(tmp_path)

        # Copy over test data
        shutil.copytree("tests/data/bids_t1w", f"{self.tmpdir}/data")
        assert filecmp.dircmp("tests/data/bids_t1w", f"{self.tmpdir}/data")

        # Create config
        self.bids_dir: str = f"{self.tmpdir}/data"
        self.pybidsdb_dir = ""
        self.pybidsdb_reset = False

    def test_pybidsdb_dir_blank(self):
        # Test non-saving (check db does not exist)
        _gen_bids_layout(
            bids_dir=self.bids_dir,
            derivatives=False,
            pybidsdb_dir=self.pybidsdb_dir,
            pybidsdb_reset=self.pybidsdb_reset,
        )
        assert not os.path.exists(self.pybidsdb_dir)

    def test_pybidsdb_dir_relative(self):
        # Update config
        self.pybidsdb_dir = "./.db"

        # Check to make sure db exists (relative path)
        _gen_bids_layout(
            bids_dir=self.bids_dir,
            derivatives=False,
            pybidsdb_dir=self.pybidsdb_dir,
            pybidsdb_reset=self.pybidsdb_reset,
        )
        assert not os.path.exists(f"{self.tmpdir}/data/.db/")

    def test_pybidsdb_dir_absolute(self):
        # Update config
        self.pybidsdb_dir = f"{self.tmpdir}/data/.db/"
        self.pybidsdb_reset = False

        # Check to make sure db exists (absolute path)
        _gen_bids_layout(
            bids_dir=self.bids_dir,
            derivatives=False,
            pybidsdb_dir=self.pybidsdb_dir,
            pybidsdb_reset=self.pybidsdb_reset,
        )
        assert os.path.exists(f"{self.tmpdir}/data/.db/")

        # Test reading of old layout when changes occur
        os.makedirs(f"{self.tmpdir}/data/sub-003/anat")
        shutil.copy(
            f"{self.bids_dir}/sub-001/anat/sub-001_acq-mprage_T1w.nii.gz",
            f"{self.bids_dir}/sub-003/anat/sub-003_acq-mprage_T1w.nii.gz",
        )
        # Check to make sure new subject not cached in layout
        layout = _gen_bids_layout(
            bids_dir=self.bids_dir,
            derivatives=False,
            pybidsdb_dir=self.pybidsdb_dir,
            pybidsdb_reset=self.pybidsdb_reset,
        )
        assert not layout.get(subject="003")

        # Test updating of layout
        self.pybidsdb_reset = True
        # Check to see if new subject in updated layout
        layout = _gen_bids_layout(
            bids_dir=self.bids_dir,
            derivatives=False,
            pybidsdb_dir=self.pybidsdb_dir,
            pybidsdb_reset=self.pybidsdb_reset,
        )
        assert layout.get(subject="003")
