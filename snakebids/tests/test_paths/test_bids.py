from __future__ import annotations

import itertools as it
import os
from pathlib import Path
from typing import Iterable

import more_itertools as itx
import pytest
from hypothesis import assume, example, given
from hypothesis import strategies as st
from pathvalidate import Platform, is_valid_filename, is_valid_filepath

from snakebids.paths import specs
from snakebids.paths.factory import bids_factory
from snakebids.paths.utils import BidsPathSpec
from snakebids.tests import strategies as sb_st
from snakebids.tests.helpers import Benchmark, is_strictly_increasing
from snakebids.utils.utils import BidsEntity


def _get_entity_tags(entities: Iterable[str]):
    return [BidsEntity.normalize(entity).tag for entity in entities]


def _values() -> st.SearchStrategy[str]:
    return st.text(min_size=1).filter(
        lambda s: is_valid_filename(s, Platform.LINUX) and "_" not in s and "-" not in s
    )


def _roots():
    return st.text().filter(lambda s: is_valid_filepath(s, Platform.LINUX))


def make_bids_testsuite(spec: BidsPathSpec):
    std_entities = {e["entity"] for e in spec}
    has_tag = {e["entity"] for e in spec if e.get("tag")}
    has_dir = {e["entity"] for e in spec if e.get("dir")}
    bids = bids_factory(spec)

    def _bids_args(
        entities: set[str] | None = std_entities,
        nonstandard: bool = True,
        custom: bool = True,
    ):
        return (
            st.dictionaries(
                keys=st.one_of(
                    # standard
                    sb_st.bids_entity(whitelist_entities=entities)
                    if entities is not None
                    else sb_st.nothing(),
                    # custom entities
                    _values()
                    .map(BidsEntity)
                    .filter(
                        lambda s: str(s)
                        not in {"datatype", "suffix", "extension", "prefix"}
                    )
                    if custom
                    else sb_st.nothing(),
                    # nonstandard entities
                    sb_st.bids_entity(
                        whitelist_entities=["datatype", "suffix", "extension"]
                    )
                    if nonstandard
                    else sb_st.nothing(),
                ),
                # The boolean here is to decide whether to use the entity or the tag in
                # the BidsEntity generated above
                values=st.tuples(
                    st.booleans(),
                    _values(),
                ),
                min_size=1,
            )
            .map(
                # we only want to use the full entity name if the spec in use
                # understands the full entity name
                lambda d: {
                    key.entity if val[0] and str(key) in has_tag else key.tag: val[1]
                    for key, val in d.items()
                }
            )
            .filter(lambda s: set(s) - {"datatype", "extension"})
        )

    class BidsTests:
        @given(_bids_args(nonstandard=False))
        def test_number_of_underscore_corresponds_to_number_entities(
            self, entities: dict[str, str]
        ):
            assert bids(**entities).count("_") == len(entities) - 1

        @given(_bids_args(nonstandard=False))
        def test_number_of_dashes_corresponds_to_number_entities(
            self, entities: dict[str, str]
        ):
            assert Path(bids(**entities)).name.count("-") == len(entities)

        @given(entities=_bids_args(), prefix=_values())
        def test_beginning_of_name_always_prefix(
            self, entities: dict[str, str], prefix: str
        ):
            assert Path(bids(prefix=prefix, **entities)).name.startswith(prefix)

        @given(
            entities=_bids_args(nonstandard=False),
            suffix=_values(),
            extension=_values(),
        )
        def test_end_of_path_always_suffix_extension(
            self, entities: dict[str, str], suffix: str, extension: str
        ):
            assert bids(suffix=suffix, extension=extension, **entities).endswith(
                suffix + extension
            )

        @given(entities=_bids_args(nonstandard=False), suffix=_values())
        def test_underscore_precedes_suffix(
            self, entities: dict[str, str], suffix: str
        ):
            assert bids(suffix=suffix, **entities)[len(suffix) * -1 - 1] == "_"

        @given(
            entities=_bids_args(nonstandard=False),
            suffix=_values() | st.none(),
            extension=_values(),
        )
        def test_underscore_does_not_precede_extension(
            self, entities: dict[str, str], suffix: str | None, extension: str
        ):
            assert (
                bids(suffix=suffix, extension=extension, **entities)[
                    len(extension) * -1 - 1
                ]
                != "_"
            )

        @given(entities=_bids_args(nonstandard=False))
        def test_no_underscore_at_end_if_no_suffix(self, entities: dict[str, str]):
            assert bids(**entities)[-1] != "_"

        @given(entities=_bids_args(), root=_roots())
        def test_beginning_of_path_always_root(
            self, entities: dict[str, str], root: str
        ):
            path = bids(root=root, **entities)
            assert path.startswith(root)
            length = len(root) - 1 if root[-1] == os.path.sep else len(root)
            assert path[length] == os.path.sep

        @example(entities={"sub": "0"}, datatype=".", root="0")
        @given(
            entities=_bids_args(nonstandard=False), datatype=_values(), root=_roots()
        )
        def test_bottom_directory_always_datatype(
            self, entities: dict[str, str], datatype: str, root: str
        ):
            # use os.path functions so that datatype=="." is treated safely
            assert datatype == os.path.basename(
                os.path.dirname(bids(root=root, datatype=datatype, **entities))
            )

        @given(entities=_bids_args(nonstandard=False), datatype=_values())
        def test_datatype_not_in_path_name(
            self, entities: dict[str, str], datatype: str
        ):
            assume(datatype not in "".join(it.chain.from_iterable(entities.items())))
            assert datatype not in Path(bids(datatype=datatype, **entities)).name

        @given(entities=_bids_args(nonstandard=False))
        def test_entities_all_in_path_as_tags(self, entities: dict[str, str]):
            tags = _get_entity_tags(entities)
            path = "_" + Path(bids(**entities)).name
            for tag in tags:
                assert f"_{tag}-" in path

        @given(entities=_bids_args(entities=has_tag, nonstandard=False, custom=False))
        def test_full_entity_names_not_in_path(self, entities: dict[str, str]):
            non_tags = [
                normed.entity
                for entity in entities
                if entity != (normed := BidsEntity.normalize(entity)).tag
            ]
            path = bids(**entities)
            for entity in non_tags:
                assert f"{entity}-" not in path

        @given(entities=_bids_args(entities=has_tag, nonstandard=False, custom=False))
        def test_long_and_short_names_cannot_be_used_simultaneously(
            self, entities: dict[str, str]
        ):
            entities = {BidsEntity.normalize(e).entity: v for e, v in entities.items()}
            tags = {BidsEntity.normalize(e).tag: v for e, v in entities.items()}
            with pytest.raises(
                ValueError,
                match="Long and short names of an entity cannot be used in the same",
            ) as err:
                bids(**entities, **tags)
            assert itx.first(entities) in err.value.args[0]
            assert itx.first(tags) in err.value.args[0]

        @given(
            entities=_bids_args(nonstandard=False, custom=False),
            custom=_bids_args(entities=None, nonstandard=False),
        )
        def test_entities_found_in_name_in_correct_order(
            self, entities: dict[str, str], custom: dict[str, str]
        ):
            # Make sure custom tags don't overlap with main tags
            assume(not set(entities) & set(custom))
            tags = _get_entity_tags(entities)
            order: list[str] = []
            for entity in spec:
                if (tag := entity.get("tag", entity["entity"])) in tags:
                    order.append(tag)
            # Custom tags come after defined tags
            order.extend(custom)
            path = "_" + Path(bids(**custom, **entities)).name
            assert set(tags) | set(custom) == set(order)
            assert is_strictly_increasing(path.index(f"_{e}-") for e in order)

        @given(
            entities=_bids_args(entities=std_entities - has_dir, nonstandard=False),
            root=_roots(),
        )
        def test_nondir_entities_dont_have_dirs(
            self, entities: dict[str, str], root: str
        ):
            assert Path(bids(root=root, **entities)).parent == Path(root)

        @given(entities=_bids_args(entities=has_dir, nonstandard=False, custom=False))
        def test_dir_entities_each_own_dir(self, entities: dict[str, str]):
            for par in itx.islice_extended(Path(bids(**entities)).parents, 0, -1):
                count = 0
                for e in list(entities):
                    tag = BidsEntity.normalize(e).tag
                    if (
                        par.name[: len(tag)] == tag
                        # if found, appears nowhere else
                        and f"{tag}-" not in str(par.parent)
                    ):
                        del entities[e]
                        count += 1
                assert count == 1
            # Check that all entities have been removed (ie found)
            assert not entities

        @given(
            entities=_bids_args(entities=has_dir, nonstandard=False, custom=False),
            root=_roots().filter(lambda s: s != "."),
        )
        def test_directories_in_correct_order(
            self, entities: dict[str, str], root: str
        ):
            tags = _get_entity_tags(entities)
            order: list[str] = []
            for entity in spec:
                if (tag := entity.get("tag", entity["entity"])) in tags:
                    order.append(tag)
            path = str(Path(bids(root=root, **entities)).parent)
            assert set(tags) == set(order)
            assert is_strictly_increasing(
                path.index(f"{os.path.sep}{e}-") for e in order
            )

        @given(entities=_bids_args(nonstandard=False), root=_roots())
        def test_values_paired_with_entities(self, entities: dict[str, str], root: str):
            path = Path(bids(root=root, **entities))
            name = "_" + path.name
            parent = os.path.sep + str(path.parent)

            def assert_follows(string: str, first: str, second: str):
                start = string.index(first) + len(first)
                assert string[start : start + len(second)] == second

            for entity, value in entities.items():
                tag = BidsEntity.normalize(entity).tag
                assert_follows(name, f"_{tag}-", value)
                if entity in has_dir:
                    assert_follows(parent, f"{os.path.sep}{tag}-", value)

        def test_bids_with_no_args_gives_empty_path(self):
            assert not bids()

        @given(
            args=st.dictionaries(
                st.sampled_from(["datatype", "prefix"]), _values(), min_size=1
            )
        )
        def test_missing_essential_entities_gives_error(self, args: dict[str, str]):
            with pytest.raises(
                ValueError,
                match="At least one of suffix, extension, or an entity must be",
            ):
                bids(**args)

    return BidsTests


TestV0_0_0 = make_bids_testsuite(specs.v0_0_0())

TestV0_10_1 = make_bids_testsuite(specs.v0_10_1())


def test_benchmark_bids(benchmark: Benchmark):
    from snakebids.paths import bids

    """If refactoring bids, be sure this benchmark doesn't needlessly increase"""
    benchmark(
        bids,
        root="foo/bar",
        subject_dir=False,
        subject="001",
        session="32",
        run="stop",
        foo="bar",
        england="britain",
        space="cosmos",
        rome="fell",
    )
