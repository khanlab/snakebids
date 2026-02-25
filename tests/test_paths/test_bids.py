from __future__ import annotations

import itertools as it
import os
import re
from collections.abc import Iterable
from pathlib import Path

import more_itertools as itx
import pytest
from hypothesis import assume, example, given
from hypothesis import strategies as st
from pathvalidate import Platform, is_valid_filename, is_valid_filepath

from snakebids.paths import specs
from snakebids.paths._factory import bids_factory
from snakebids.paths._utils import OPTIONAL_WILDCARD, BidsPathSpec
from snakebids.snakemake_compat import regex_from_filepattern
from snakebids.utils.snakemake_templates import SnakemakeFormatter
from snakebids.utils.utils import BidsEntity
from tests import strategies as sb_st
from tests.helpers import Benchmark, is_strictly_increasing
from tests.test_snakemake_templates.strategies import safe_field_names


def _get_entity_tags(entities: Iterable[str]):
    return [BidsEntity.normalize(entity).tag for entity in entities]


def _values() -> st.SearchStrategy[str]:
    return st.text(min_size=1).filter(
        lambda s: is_valid_filename(s, Platform.LINUX) and not (set("_-.") & set(s))
    )


def _field_names():
    return safe_field_names(min_size=1).filter(lambda s: not s.isdigit())


def _identifiers():
    return st.from_regex(r"[a-zA-Z][a-zA-Z0-9]*", fullmatch=True)


def _roots():
    return st.text().filter(lambda s: is_valid_filepath(s, Platform.LINUX))


def make_bids_testsuite(spec: BidsPathSpec):
    std_entities = {e["entity"] for e in spec}
    has_tag = {e["entity"] for e in spec if e.get("tag")}
    has_dir = {e["entity"] for e in spec if e.get("dir")}
    bids = bids_factory(spec)

    def _entities(
        entities: set[str] | None = std_entities,
        nonstandard: bool = True,
        custom: st.SearchStrategy[str] | None = None,
        prefix: bool = False,
    ):
        std_ents = (
            sb_st.bids_entity(whitelist_entities=entities)
            if entities is not None
            else sb_st.nothing()
        )
        custom_ents = (
            (_values() if custom is None else custom)
            .map(BidsEntity.from_tag)
            .filter(
                lambda s: str(s)
                not in (
                    {"datatype", "suffix", "extension", "prefix"} | set(std_entities)
                )
            )
        )
        nonstd_ents = (
            sb_st.bids_entity(whitelist_entities=["datatype", "suffix", "extension"])
            if nonstandard
            else sb_st.nothing()
        )
        prefix_ent = st.just(BidsEntity("prefix")) if prefix else sb_st.nothing()
        return std_ents | custom_ents | nonstd_ents | prefix_ent

    def _bids_args(
        entities: set[str] | None = std_entities,
        nonstandard: bool = True,
        custom: st.SearchStrategy[str] | None = None,
        prefix: bool = False,
        use_wildcard_only: bool = False,
    ):
        def normalize(d: dict[BidsEntity, tuple[bool, str]]):
            result: dict[str, str] = {}
            for key, val in d.items():
                if use_wildcard_only:
                    newkey = key.wildcard
                elif val[0] and str(key) in has_tag:
                    newkey = key.entity
                else:
                    newkey = key.tag
                newval = val[1]
                if newkey == "extension":
                    newval = f".{newval}"

                result[newkey] = newval
            return result

        return (
            st.dictionaries(
                keys=_entities(entities, nonstandard, custom, prefix),
                # The boolean here is to decide whether to use the entity or the tag in
                # the BidsEntity generated above
                values=st.tuples(
                    st.booleans(),
                    _values(),
                ),
                min_size=1,
            )
            .map(normalize)
            .filter(lambda s: set(s) - {"datatype", "extension"})
        )

    class BidsTests:
        @given(_bids_args(nonstandard=False))
        def test_number_of_underscore_corresponds_to_number_entities(
            self, entities: dict[str, str]
        ):
            assert bids(**entities).count("_") == len(entities) - 1

        @given(_bids_args(nonstandard=False))
        def test_underscores_never_doubled(self, entities: dict[str, str]):
            assert "__" not in bids(**entities)

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

        @given(entities=_bids_args())
        def test_never_underscore_at_name_beginning(self, entities: dict[str, str]):
            assert Path(bids(**entities)).name[0] != "_"

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

        @given(
            entities=_bids_args(
                entities=has_tag, nonstandard=False, custom=st.nothing()
            )
        )
        def test_full_entity_names_not_in_path(self, entities: dict[str, str]):
            non_tags = [
                normed.entity
                for entity in entities
                if entity != (normed := BidsEntity.normalize(entity)).tag
            ]
            path = bids(**entities)
            for entity in non_tags:
                assert f"{entity}-" not in path

        @given(
            entities=_bids_args(
                entities=has_tag, nonstandard=False, custom=st.nothing()
            )
        )
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
            entities=_bids_args(nonstandard=False, custom=st.nothing()),
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

        @given(
            entities=_bids_args(
                entities=has_dir, nonstandard=False, custom=st.nothing()
            )
        )
        def test_dir_entities_each_own_dir(self, entities: dict[str, str]):
            for par in Path(bids(**entities)).parents[:-1]:
                count = 0
                for e in list(entities):
                    tag = BidsEntity.normalize(e).tag
                    if (
                        # tag found in parent
                        par.name[: len(tag)] == tag
                        # and nowhere else
                        and f"{tag}-" not in str(par.parent)
                    ):
                        del entities[e]
                        count += 1
                assert count == 1
            # Check that all entities have been removed (ie found)
            assert not entities

        @given(
            entities=_bids_args(
                entities=has_dir, nonstandard=False, custom=st.nothing()
            ),
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

        @given(entities=_bids_args(nonstandard=False))
        def test_values_paired_with_entities(self, entities: dict[str, str]):
            path = Path(bids(**entities))

            name_mapping = dict(s.split("-") for s in path.name.split("_"))
            parent_mapping = dict(s.split("-") for s in path.parent.parts)

            for entity, value in entities.items():
                tag = BidsEntity.normalize(entity).tag
                assert name_mapping[tag] == value
                if entity in has_dir:
                    assert parent_mapping[tag] == value

        def test_bids_with_no_args_gives_empty_path(self):
            assert bids() == ""

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

        @given(entities=_bids_args(nonstandard=False), suffix=st.text())
        def test_suffix_at_end(self, entities: dict[str, str], suffix: str):
            assert bids(suffix=suffix, **entities).endswith(suffix)

        @given(
            entities=_bids_args(nonstandard=False),
            suffix=st.text(),
            ext=st.text(),
        )
        def test_extension_at_end(
            self, entities: dict[str, str], suffix: str, ext: str
        ):
            assert bids(suffix=suffix, extension=ext, **entities).endswith(ext)

        @given(
            mandatory=_bids_args(custom=_field_names(), use_wildcard_only=True),
            optional=_bids_args(custom=_field_names(), use_wildcard_only=True),
        )
        def test_optional_wildcards_format_into_bids_paths(
            self, mandatory: dict[str, str], optional: dict[str, str]
        ):
            """Test that optional wildcards compose correctly with hypothesis."""
            args = mandatory | optional
            reference = bids(**args)
            template = bids(
                root="",
                **(
                    {
                        k: v.replace("{", "{{").replace("}", "}}")
                        for k, v in mandatory.items()
                    }
                    | dict.fromkeys(optional, OPTIONAL_WILDCARD)
                ),
            )

            formatter = SnakemakeFormatter()

            assert formatter.format(template, **args) == reference

        @given(
            mandatory=_bids_args(
                use_wildcard_only=True, custom=_identifiers(), prefix=True
            ),
            optional=_bids_args(use_wildcard_only=True, custom=_identifiers()),
        )
        def test_optional_wildcards_match_bids_paths(
            self, mandatory: dict[str, str], optional: dict[str, str]
        ):
            """Test that optional wildcards compose correctly with hypothesis."""
            mandatory = {k: v.replace("{", "{{") for k, v in mandatory.items()}
            args = mandatory | optional
            reference = bids(**args)
            template = bids(
                root="",
                **(mandatory | dict.fromkeys(optional, OPTIONAL_WILDCARD)),
            )
            regex = regex_from_filepattern(template)
            match = re.match(regex, reference)
            assert match is not None
            gathered_args = match.groupdict()
            for arg, val in optional.items():
                assert gathered_args[arg] == val

        @given(
            mandatory=_bids_args(use_wildcard_only=True, custom=_identifiers()),
            optional=_bids_args(use_wildcard_only=True, custom=_identifiers()),
        )
        def test_pathlib_roundtrip(
            self, mandatory: dict[str, str], optional: dict[str, str]
        ):
            template = bids(
                root="",
                **(mandatory | dict.fromkeys(optional, OPTIONAL_WILDCARD)),
            )
            assert str(Path(template)) == template

        def test_prefix_may_not_be_optional(self):
            with pytest.raises(
                ValueError,
                match="prefix may not be specified as optional",
            ):
                bids(prefix=OPTIONAL_WILDCARD)

    return BidsTests


TestV0_0_0 = make_bids_testsuite(specs.v0_0_0())

TestV0_11_0 = make_bids_testsuite(specs.v0_11_0())

TestV0_15_0 = make_bids_testsuite(specs.v0_15_0())


class TestBenchmarkBids:
    call = {  # noqa: RUF012
        "root": "foo/bar",
        "subject_dir": False,
        "subject": "001",
        "session": "32",
        "run": "stop",
        "foo": "bar",
        "england": "britain",
        "space": "cosmos",
        "rome": "fell",
        "suffix": "suffix",
        "extension": ".ext",
    }

    def test_new_bids(self, benchmark: Benchmark):
        """If refactoring bids, be sure this benchmark doesn't needlessly increase"""
        new_bids = bids_factory(specs.v0_0_0())
        benchmark(new_bids, **self.call)  # type: ignore
