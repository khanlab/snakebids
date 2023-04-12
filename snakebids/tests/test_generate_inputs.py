import filecmp
import functools as ft
import itertools as it
import keyword
import os
import re
import shutil
import tempfile
from collections import defaultdict
from copy import deepcopy
from pathlib import Path
from typing import Dict, Iterable, List, NamedTuple, Optional, TypeVar

import attrs
import more_itertools as itx
import pytest
from bids import BIDSLayout
from hypothesis import HealthCheck, assume, example, given, settings
from hypothesis import strategies as st
from pyfakefs.fake_filesystem import FakeFilesystem
from pytest_mock import MockerFixture

from snakebids.core.construct_bids import bids
from snakebids.core.datasets import BidsComponent, BidsDataset
from snakebids.core.input_generation import (
    _all_custom_paths,
    _gen_bids_layout,
    _generate_filters,
    _get_lists_from_bids,
    _parse_custom_path,
    generate_inputs,
)
from snakebids.exceptions import ConfigError, PybidsError
from snakebids.tests import strategies as sb_st
from snakebids.tests.helpers import (
    BidsListCompare,
    allow_tmpdir,
    create_dataset,
    debug,
    get_bids_path,
    get_zip_list,
    reindex_dataset,
)
from snakebids.types import InputsConfig
from snakebids.utils import sb_itertools as sb_it
from snakebids.utils.utils import BidsEntity

T = TypeVar("T")


class TestFilterBools:
    @pytest.fixture(autouse=True)
    def bids_fs(self, bids_fs: Optional[FakeFilesystem]):
        return bids_fs

    @pytest.fixture
    def tmpdir(self, fakefs_tmpdir: Path):
        return fakefs_tmpdir

    def disambiguate_components(self, dataset: BidsDataset):
        assert len(dataset) == 2
        comp1, comp2 = dataset.values()
        return sorted([comp1, comp2], key=lambda comp: len(comp.entities.keys()))

    def get_extra_entity(self, dataset: BidsDataset) -> str:
        return itx.one(
            ft.reduce(
                lambda set_, comp: set_.symmetric_difference(comp.entities.keys()),
                dataset.values(),
                set(),
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
                    for wildcard in shorter.input_wildcards.keys()
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
                    for wildcard in longer.input_wildcards.keys()
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
                    for wildcard in shorter.input_wildcards.keys()
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
    @given(dataset=sb_st.datasets())
    def test_entity_excluded_when_filter_true(self, tmpdir: Path, dataset: BidsDataset):
        root = tempfile.mkdtemp(dir=tmpdir)
        create_dataset(root, dataset)
        _, longer = self.disambiguate_components(dataset)
        extra_entity = self.get_extra_entity(dataset)
        pybids_inputs: InputsConfig = {
            longer.input_name: {
                "wildcards": [
                    BidsEntity.from_tag(wildcard).entity
                    for wildcard in longer.input_wildcards.keys()
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


class TestAbsentConfigEntries:
    def get_entities(self, root):
        # Generate directory
        entities = {"subject": ["001", "002"], "acq": sorted(["foo", "bar"])}
        zip_list = defaultdict(list)
        for e in it.product(*entities.values()):
            d = dict(zip(entities.keys(), e))
            for key, val in d.items():
                zip_list[key].append(val)
            path = Path(bids(root, datatype="anat", suffix="T1w.nii.gz", **d))
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()
        return entities, zip_list

    def test_missing_filters(self, tmpdir: Path):
        entities, zip_list = self.get_entities(tmpdir)

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
        template = BidsDataset({"t1": BidsComponent("t1", config["t1"].path, zip_list)})
        # Order of the subjects is not deterministic
        assert template == config
        assert config.subj_wildcards == {"subject": "{subject}"}

    def test_missing_wildcards(self, tmpdir: Path):
        entities, zip_list = self.get_entities(tmpdir)

        # create config
        derivatives = False
        pybids_inputs: InputsConfig = {
            "t1": {
                "filters": {"acquisition": "foo", "subject": "001"},
            }
        }

        # Simplest case -- one input type, using pybids
        config = generate_inputs(
            pybids_inputs=pybids_inputs,
            bids_dir=tmpdir,
            derivatives=derivatives,
            pybids_config=str(Path(__file__).parent / "data" / "custom_config.json"),
        )
        template = BidsDataset({"t1": BidsComponent("t1", config["t1"].path, {})})
        assert template == config
        assert config.subj_wildcards == {"subject": "{subject}"}


class TestGenerateFilter:
    valid_chars = st.characters(blacklist_characters=["\n"])
    st_lists_or_text = st.lists(st.text(valid_chars)) | st.text(valid_chars)

    @given(st.tuples(st_lists_or_text, st_lists_or_text))
    def test_throws_error_if_labels_and_excludes_are_given(self, args):
        with pytest.raises(ValueError):
            _generate_filters(*args)

    @given(st_lists_or_text)
    def test_returns_participant_label_as_list(self, label):
        result = _generate_filters(label)[0]
        if isinstance(label, str):
            assert result == [label]
        else:
            assert result == label

    @given(
        st_lists_or_text,
        st.lists(st.text(valid_chars, min_size=1), min_size=1),
        st.text(valid_chars, min_size=1, max_size=3),
    )
    def test_exclude_gives_regex_that_matches_anything_except_exclude(
        self, excluded, dummy_values, padding
    ):
        # Make sure the dummy_values and padding we'll be testing against are different
        # from our test values
        for value in dummy_values:
            assume(value not in itx.always_iterable(excluded))
        assume(padding not in itx.always_iterable(excluded))

        result = _generate_filters(exclude=excluded)
        assert result[1] is True
        assert isinstance(result[0], list)
        assert len(result[0]) == 1

        # We match any value that isn't the exclude string
        for value in dummy_values:
            assert re.match(result[0][0], value)

        for exclude in itx.always_iterable(excluded):
            # We don't match the exclude string
            assert re.match(result[0][0], exclude) is None

            # Addition of random strings before and/or after lets the match occur again
            assert re.match(result[0][0], padding + exclude)
            assert re.match(result[0][0], exclude + padding)
            assert re.match(result[0][0], padding + exclude + padding)


PathEntities = NamedTuple(
    "PathEntities",
    [
        ("entities", Dict[str, List[str]]),
        ("template", Path),
        ("filters", Dict[str, List[str]]),
    ],
)


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
    # We need to explicitely exclude keywords here because the current implementation
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

    def get_subset(of: Iterable[T]) -> List[T]:
        return draw(
            st.lists(st.sampled_from([*of]), unique=True, max_size=len(entities))
        )

    # Collect dir_entities and filters
    dir_entities = get_subset(entities.keys())
    filtered_entities = get_subset(entities.keys())
    filter_selections = [get_subset(entities[entity]) for entity in filtered_entities]
    filters = dict(zip(filtered_entities, filter_selections))

    # Compose the path template
    dir_template = Path(*(f"{{{entity}}}" for entity in dir_entities))
    name_template = "_".join(f"{entity}-{{{entity}}}" for entity in entities)
    template = dir_template / name_template

    return PathEntities(entities, template, filters)


class TestCustomPaths:
    @pytest.fixture()
    def temp_dir(self, fakefs_tmpdir: Path):
        return fakefs_tmpdir

    def generate_test_directory(
        self, entities: Dict[str, List[str]], template: Path, tmpdir: Path
    ):
        root = Path(tempfile.mkdtemp(prefix="hypothesis-", dir=tmpdir))
        # Generate fake directory structure
        for values in it.product(*entities.values()):
            name_value = dict(zip(entities.keys(), values))
            path = str(template).format(**name_value)
            (root / path).parent.mkdir(parents=True, exist_ok=True)
            (root / path).touch()
        return root / template

    def test_benchmark_test_custom_paths(self, benchmark, tmp_path: Path):
        entities = {"A": ["A", "B", "C"], "B": ["1", "2", "3"]}
        template = Path("{A}/A-{A}_B-{B}")
        test_path = self.generate_test_directory(entities, template, tmp_path)
        benchmark(_parse_custom_path, test_path)

    @allow_tmpdir
    @given(path_entities=path_entities())
    def test_collects_all_paths_when_no_filters(
        self,
        path_entities: PathEntities,
        temp_dir: Path,
    ):
        entities, template, _ = path_entities
        test_path = self.generate_test_directory(entities, template, temp_dir)

        # Test without any filters
        result = _parse_custom_path(test_path)
        zip_lists = get_zip_list(entities, it.product(*entities.values()))
        assert BidsComponent(
            "foo", get_bids_path(zip_lists), zip_lists
        ) == BidsComponent("foo", get_bids_path(result), result)

    @allow_tmpdir
    @given(path_entities=path_entities())
    def test_collects_only_filtered_entities(
        self,
        path_entities: PathEntities,
        temp_dir: Path,
    ):
        entities, template, filters = path_entities
        test_path = self.generate_test_directory(entities, template, temp_dir)

        # Test with filters
        result_filtered = _parse_custom_path(test_path, **filters)
        zip_lists = {
            # Start with empty lists for each key, otherwise keys will be missing
            **{key: [] for key in entities},
            # Override entities with relevant filters before making zip lists
            **get_zip_list(entities, it.product(*{**entities, **filters}.values())),
        }
        assert BidsComponent(
            "foo", get_bids_path(zip_lists), zip_lists
        ) == BidsComponent("foo", get_bids_path(result_filtered), result_filtered)

    @allow_tmpdir
    @given(path_entities=path_entities())
    def test_collect_all_but_filters_when_exclusion_filters_used(
        self,
        path_entities: PathEntities,
        temp_dir: Path,
    ):
        entities, template, filters = path_entities
        test_path = self.generate_test_directory(entities, template, temp_dir)
        # Test with exclusion filters
        exclude_filters = {
            # We use _generate_filter to get our exclusion regex. This function was
            # tested previously
            key: _generate_filters(exclude=values)[0]
            for key, values in filters.items()
        }
        result_excluded = _parse_custom_path(
            test_path, regex_search=True, **exclude_filters
        )

        entities_excluded = {
            entity: [value for value in values if value not in filters.get(entity, [])]
            for entity, values in entities.items()
        }
        zip_lists = {
            # Start with empty lists for each key, otherwise keys will be missing
            **{key: [] for key in entities},
            # Override entities with relevant filters before making zip lists
            **get_zip_list(entities, it.product(*entities_excluded.values())),
        }

        assert BidsComponent(
            "foo", get_bids_path(zip_lists), zip_lists
        ) == BidsComponent("foo", get_bids_path(result_excluded), result_excluded)


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
        pybids_config=Path(__file__).parent / "data" / "custom_config.json",
    )
    template = BidsDataset(
        {
            "t1": BidsComponent(
                "t1",
                bids(
                    tmpdir,
                    datatype="anat",
                    subject="{subject}",
                    foo="{foo}",
                    suffix="T1w.nii.gz",
                ),
                {"foo": ["0", "1"], "subject": ["001", "001"]},
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
                Path(__file__).parent / "data" / "custom_config_nonstandard.json"
            ),
        )


def test_t1w():
    # create config
    real_bids_dir = "snakebids/tests/data/bids_t1w"
    derivatives = False
    pybids_inputs: InputsConfig = {
        "t1": {
            "filters": {"suffix": "T1w"},
            "wildcards": ["acquisition", "subject", "session", "run"],
        }
    }

    # Can't define particpant_label and exclude_participant_label
    with pytest.raises(ValueError) as v_error:
        result = generate_inputs(
            pybids_inputs=pybids_inputs,
            bids_dir=real_bids_dir,
            derivatives=derivatives,
            participant_label="001",
            exclude_participant_label="002",
        )
    assert v_error.value.args[0] == (
        "Cannot define both participant_label and "
        "exclude_participant_label at the same time"
    )

    # Simplest case -- one input type, using pybids
    result = generate_inputs(
        pybids_inputs=pybids_inputs,
        bids_dir=real_bids_dir,
        derivatives=derivatives,
    )
    template = BidsDataset(
        {
            "t1": BidsComponent(
                "t1",
                result["t1"].path,
                {"acq": ["mprage", "mprage"], "subject": ["001", "002"]},
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
                "scan",
                result["scan"].path,
                {
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
                    "t1",
                    result["t1"].path,
                    {
                        "acq": ["mprage", "mprage"],
                        "subject": ["001", "002"],
                    },
                ),
                "t2": BidsComponent("t2", result["t2"].path, {"subject": ["002"]}),
            }
        )
        # Order of the subjects is not deterministic
        assert result.subjects in [["001", "002"], ["002", "001"]]

        assert result.sessions == []
        assert result.subj_wildcards == {"subject": "{subject}"}


def test_t1w_with_dict():
    # create config
    real_bids_dir = "snakebids/tests/data/bids_t1w"
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


def test_get_lists_from_bids_raises_pybids_error():
    """Test that we wrap a cryptic AttributeError from pybids with PybidsError.

    Pybids raises an AttributeError when a BIDSLayout.get is called with scope not
    equal to 'all' on a layout that indexes a dataset without a
    dataset_description.json. We wrap this error with something a bit less cryptic, so
    this test ensures that that behaviour is still present.
    """
    layout = BIDSLayout("snakebids/tests/data/bids_t1w", validate=False)
    with pytest.raises(PybidsError):
        next(_get_lists_from_bids(layout, {"t1": {"filters": {"scope": "raw"}}}))


def test_get_lists_from_bids():
    bids_dir = "snakebids/tests/data/bids_t1w"
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

        result = _get_lists_from_bids(layout, pybids_inputs)
        for bids_lists in result:
            if bids_lists.input_name == "t1":
                template = BidsComponent(
                    "t1",
                    wildcard_path_t1,
                    {
                        "acq": ["mprage", "mprage"],
                        "subject": ["001", "002"],
                    },
                )
                assert template == bids_lists
            elif bids_lists.input_name == "t2":
                assert bids_lists.input_path == wildcard_path_t2
                template = BidsComponent(
                    "t2",
                    wildcard_path_t2,
                    {
                        "subject": ["002"],
                    },
                )
                assert template == bids_lists


class TestGenBidsLayout:
    @pytest.fixture
    def tmpdir(self, fakefs_tmpdir: Path):
        return fakefs_tmpdir

    @pytest.fixture(autouse=True)
    def bids_fs(self, bids_fs: Optional[FakeFilesystem]):
        return bids_fs

    def test_gen_layout_returns_valid_dataset(self, tmpdir: Path):
        dataset = sb_st.datasets().example()
        create_dataset(tmpdir, dataset)
        assert _gen_bids_layout(tmpdir, False, None, False, None)

    def test_invalid_path_raises_error(self, tmpdir: Path):
        with pytest.raises(ValueError):
            _gen_bids_layout(tmpdir / "foo", False, None, False)


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


@settings(
    deadline=800,
    suppress_health_check=[
        HealthCheck.function_scoped_fixture,
        HealthCheck.too_slow,
    ],
)
@given(data=st.data())
def test_generate_inputs(data: st.DataObject, bids_fs: Path, fakefs_tmpdir: Path):
    root = tempfile.mkdtemp(dir=fakefs_tmpdir)
    dataset = data.draw(sb_st.datasets(root=Path(root)))
    reindexed = reindex_dataset(root, dataset)
    assert reindexed == dataset
    assert reindexed.layout is not None


def test_when_all_custom_paths_no_layout_indexed(
    bids_fs: Path, fakefs_tmpdir: Path, mocker: MockerFixture
):
    root = tempfile.mkdtemp(dir=fakefs_tmpdir)
    # The content of the dataset is irrelevant to this test, so one example suffices
    dataset = sb_st.datasets_one_comp(root=Path(root)).example()

    spy = mocker.spy(BIDSLayout, "__init__")
    reindexed = reindex_dataset(root, dataset, use_custom_paths=True)
    assert reindexed == dataset
    assert reindexed.layout is None
    spy.assert_not_called()


class TestDB:
    @pytest.fixture(autouse=True)
    def start(self, tmpdir):
        self.tmpdir = tmpdir.strpath

        # Copy over test data
        shutil.copytree("snakebids/tests/data/bids_t1w", f"{self.tmpdir}/data")
        assert filecmp.dircmp("snakebids/tests/data/bids_t1w", f"{self.tmpdir}/data")

        # Create config
        self.bids_dir = f"{self.tmpdir}/data"
        self.pybids_db = {"database_dir": "", "reset_database": False}

    def test_database_dir_blank(self):
        # Test non-saving (check db does not exist)
        _gen_bids_layout(
            bids_dir=self.bids_dir,
            derivatives=False,
            pybids_database_dir=self.pybids_db.get("database_dir"),
            pybids_reset_database=self.pybids_db.get("reset_database"),
        )
        assert not os.path.exists(self.pybids_db.get("database_dir"))

    def test_database_dir_relative(self):
        # Update config
        self.pybids_db["database_dir"] = "./.db"

        # Check to make sure db exists (relative path)
        _gen_bids_layout(
            bids_dir=self.bids_dir,
            derivatives=False,
            pybids_database_dir=self.pybids_db.get("database_dir"),
            pybids_reset_database=self.pybids_db.get("reset_database"),
        )
        assert not os.path.exists(f"{self.tmpdir}/data/.db/")

    def test_database_dir_absolute(self):
        # Update config
        self.pybids_db["database_dir"] = f"{self.tmpdir}/data/.db/"
        self.pybids_db["reset_database"] = False

        # Check to make sure db exists (absolute path)
        _gen_bids_layout(
            bids_dir=self.bids_dir,
            derivatives=False,
            pybids_database_dir=self.pybids_db.get("database_dir"),
            pybids_reset_database=self.pybids_db.get("reset_database"),
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
            pybids_database_dir=self.pybids_db.get("database_dir"),
            pybids_reset_database=self.pybids_db.get("reset_database"),
        )
        assert not layout.get(subject="003")

        # Test updating of layout
        self.pybids_db["reset_database"] = True
        # Check to see if new subject in updated layout
        layout = _gen_bids_layout(
            bids_dir=self.bids_dir,
            derivatives=False,
            pybids_database_dir=self.pybids_db.get("database_dir"),
            pybids_reset_database=self.pybids_db.get("reset_database"),
        )
        assert layout.get(subject="003")
