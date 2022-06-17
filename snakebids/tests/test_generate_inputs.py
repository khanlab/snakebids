from __future__ import absolute_import

import copy
import filecmp
import itertools as it
import keyword
import os
import re
import shutil
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, NamedTuple, TypeVar

import more_itertools as itx
import pytest
from bids import BIDSLayout
from hypothesis import assume, given
from hypothesis import strategies as st

from snakebids.core.construct_bids import bids
from snakebids.core.input_generation import (
    BidsComponent,
    BidsDataset,
    _gen_bids_layout,
    _generate_filters,
    _get_lists_from_bids,
    _parse_custom_path,
    generate_inputs,
)
from snakebids.tests import strategies as sb_st
from snakebids.tests.helpers import BidsListCompare, get_zip_list, setify

T = TypeVar("T")


class TestBidsInputEq:
    @given(sb_st.bids_input(), sb_st.everything_except(BidsComponent))
    def test_other_types_are_unequal(self, input: BidsComponent, other: Any):
        assert input != other

    def test_empty_BidsInput_are_equal(self):
        assert BidsComponent("", "", {}) == BidsComponent("", "", {})
        assert BidsComponent("", "", {"foo": [], "bar": []}) == BidsComponent(
            "", "", {"foo": [], "bar": []}
        )

    @given(sb_st.bids_input())
    def test_copies_are_equal(self, input: BidsComponent):
        cp = copy.deepcopy(input)
        assert cp == input

    @given(sb_st.bids_input())
    def test_mutation_makes_unequal(self, input: BidsComponent):
        cp = copy.deepcopy(input)
        next(iter(cp.input_zip_lists.values()))[0] += "foo"
        assert cp != input

    @given(sb_st.bids_input(), st.data())
    def test_extra_entities_makes_unequal(
        self, input: BidsComponent, data: st.DataObject
    ):
        cp = copy.deepcopy(input)
        new_entity = data.draw(
            sb_st.bids_value().filter(lambda s: s not in input.input_zip_lists)
        )
        cp.input_zip_lists[new_entity] = []
        next(iter(cp.input_zip_lists.values()))[0] += "foo"
        assert cp != input

    @given(sb_st.bids_input())
    def test_order_doesnt_affect_equality(self, input: BidsComponent):
        cp = copy.deepcopy(input)
        for l in cp.input_zip_lists:
            cp.input_zip_lists[l].reverse()
        assert cp == input

    @given(sb_st.bids_input())
    def test_paths_must_be_identical(self, input: BidsComponent):
        cp = BidsComponent(
            input.input_name, input.input_path + "foo", input.input_zip_lists
        )
        assert cp != input


class TestBidsComponentProperties:
    @given(st.data(), st.integers(min_value=1, max_value=2))
    def test_input_lists_derives_from_zip_lists(
        self, data: st.DataObject, min_size: int
    ):
        input_lists: Dict[str, List[str]] = data.draw(
            sb_st.bids_input_lists(min_size, max_size=5)
        )

        # Due to the product, we can delete some of the combinations and still
        # regenerate our input_lists
        combs = list(it.product(*input_lists.values()))[min_size - 1 :]
        zip_lists = get_zip_list(input_lists, combs)
        path = sb_st.get_bids_path(zip_lists)

        assert setify(
            BidsComponent(
                input_name="foo", input_path=path, input_zip_lists=zip_lists
            ).input_lists
        ) == setify(input_lists)

    @given(st.dictionaries(sb_st.bids_entity(), sb_st.bids_value(), min_size=1))
    def test_input_wildcards_derives_from_zip_lists(
        self,
        bids_entities: Dict[str, str],
    ):
        zip_lists = {entity: [val] for entity, val in bids_entities.items()}
        bids_input = BidsComponent(
            input_name="foo", input_path="foo", input_zip_lists=zip_lists
        )

        wildstr = ".".join(bids_input.input_wildcards.values())
        first = wildstr.format(**bids_input.input_wildcards)
        second = first.format(**bids_entities)
        assert set(second.split(".")) == set(bids_entities.values())


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
        pybids_inputs = {
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
            use_bids_inputs=True,
        )
        template = BidsDataset(
            {"t1": BidsComponent("t1", config.input_path["t1"], zip_list)}
        )
        # Order of the subjects is not deterministic
        assert template == config
        assert config.subj_wildcards == {"subject": "{subject}"}

    def test_missing_wildcards(self, tmpdir: Path):
        entities, zip_list = self.get_entities(tmpdir)

        # create config
        derivatives = False
        pybids_inputs = {
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
            use_bids_inputs=True,
        )
        template = BidsDataset({"t1": BidsComponent("t1", config.input_path["t1"], {})})
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
    @pytest.fixture(scope="class")
    def temp_dir(self, tmp_path_factory: pytest.TempPathFactory):
        return tmp_path_factory.mktemp("test-custom-paths-")

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
        assert BidsComponent("foo", "foo", zip_lists) == BidsComponent(
            "foo", "foo", result
        )

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
        assert BidsComponent("foo", "foo", zip_lists) == BidsComponent(
            "foo", "foo", result_filtered
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

        assert BidsComponent("foo", "foo", zip_lists) == BidsComponent(
            "foo", "foo", result_excluded
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
    pybids_inputs = {
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
        pybids_config=(Path(__file__).parent / "data" / "custom_config.json"),
        use_bids_inputs=True,
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
    assert result.input_wildcards == {"t1": {"foo": "{foo}", "subject": "{subject}"}}
    # Order of the subjects is not deterministic
    assert result.subj_wildcards == {"subject": "{subject}"}


def test_t1w():
    # create config
    real_bids_dir = "snakebids/tests/data/bids_t1w"
    derivatives = False
    pybids_inputs = {
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
        use_bids_inputs=True,
    )
    template = BidsDataset(
        {
            "t1": BidsComponent(
                "t1",
                result.input_path["t1"],
                {"acq": ["mprage", "mprage"], "subject": ["001", "002"]},
            )
        }
    )
    assert template == result

    # Order of the subjects is not deterministic
    assert result.subjects in [["001", "002"], ["002", "001"]]
    assert result.sessions == []
    assert result.subj_wildcards == {"subject": "{subject}"}

    pybids_inputs_suffix = {
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
        use_bids_inputs=True,
    )
    assert result.input_lists == {
        "scan": {"acq": ["mprage"], "subject": ["001"], "suffix": ["T1w"]}
    }
    template = BidsDataset(
        {
            "scan": BidsComponent(
                "scan",
                result.input_path["scan"],
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
            use_bids_inputs=True,
        )
        template = BidsDataset(
            {
                "t1": BidsComponent(
                    "t1",
                    result.input_path["t1"],
                    {
                        "acq": ["mprage", "mprage"],
                        "subject": ["001", "002"],
                    },
                ),
                "t2": BidsComponent(
                    "t2", result.input_path["t2"], {"subject": ["002"]}
                ),
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
    pybids_inputs = {
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

    pybids_inputs_suffix = {
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
