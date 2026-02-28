"""Tests for snakemake template formatting utilities."""

from __future__ import annotations

import string

import hypothesis.strategies as st
import more_itertools as itx
import pytest
from hypothesis import assume, example, given

import tests.strategies as sb_st
from snakebids import bids
from snakebids.paths import OPTIONAL_WILDCARD
from snakebids.utils.snakemake_templates import (
    SnakemakeFormatter,
    SnakemakeWildcards,
    _HAS_RUST_PARSE,
)
from tests.helpers import Benchmark
from tests.test_snakemake_templates.strategies import (
    constraints,
    field_names,
    literals,
    safe_field_names,
)


class TestParserBenchmarks:
    times = 10000
    text = "fewaigq{{eafwefao{wafwa,waffwahaw,a:fa}{fea!f:wagwa}jawjgaw{fwa,afewaew}k"

    @staticmethod
    def run(
        formatter: string.Formatter,
        text: str,
    ):
        return list(formatter.parse(text))

    def test_benchmark_native_formatter(self, benchmark: Benchmark):
        assert benchmark(
            self.run,
            string.Formatter(),
            self.text * self.times,
        )

    def test_benchmark_custom_formatter(self, benchmark: Benchmark):
        s = self.text * self.times
        assert benchmark(self.run, SnakemakeFormatter(), s)

    def test_benchmark_rust_formatter(self, benchmark: Benchmark):
        s = self.text * self.times
        # self.run(SnakemakeFormatter(use_rust=True), s)
        assert benchmark(self.run, SnakemakeFormatter(use_rust=True), s)


class TestParse:
    """Tests for SnakemakeFormatter.parse() method."""

    @example(literals=["}}"], wildcards=["{}", "{}"])
    @given(
        literals=st.lists(literals()),
        wildcards=st.lists(
            field_names(exclude_characters="!:").map(
                lambda s: f"{{{s}}}".replace("[", "[]")
            )
        ),
    )
    def test_custom_parse_equivalent_to_native(
        self, literals: list[str], wildcards: list[str]
    ):
        path = "".join(itx.interleave_longest(literals, wildcards))
        assert list(string.Formatter().parse(path)) == list(
            SnakemakeFormatter().parse(path)
        )

    @given(name=field_names(), constraint=constraints())
    def test_parse_strips_constraints_from_field_names(
        self, name: str, constraint: str
    ):
        """Test that parse() strips constraints from field names."""
        formatter = SnakemakeFormatter()
        result = next(formatter.parse(f"{{{name},{constraint}}}"))[1]

        assert result == name

    @given(name=field_names(), format_spec=field_names())
    def test_parse_does_not_support_format_spec(self, name: str, format_spec: str):
        formatter = SnakemakeFormatter()
        template = f"{{{name}:{format_spec}}}"
        result = next(formatter.parse(template))

        # format_spec should be empty string
        assert result[1] == template[1:-1]
        assert result[2] == ""

    @given(name=field_names(), conversion=field_names())
    def test_parse_does_not_support_conversion(self, name: str, conversion: str):
        formatter = SnakemakeFormatter()
        template = f"{{{name}!{conversion}}}"
        result = next(formatter.parse(template))

        # conversion should be None
        assert result[1] == template[1:-1]
        assert result[3] is None

    def test_parse_raises_error_for_missing_closing_brace(self):
        """Test that raises ValueError for missing closing brace."""
        formatter = SnakemakeFormatter()
        with pytest.raises(ValueError, match="expected '}' before end of string"):
            list(formatter.parse("prefix_{subject"))

    def test_parse_raises_error_for_trailing_opening_brace(self):
        formatter = SnakemakeFormatter()
        with pytest.raises(ValueError, match="expected '}' before end of string"):
            list(formatter.parse("prefix_{"))

    def test_parse_raises_error_for_trailing_closing_brace(self):
        formatter = SnakemakeFormatter()
        with pytest.raises(ValueError, match="unexpected '}' in string"):
            list(formatter.parse("prefix_}"))

    def test_parse_raises_error_for_nested_opening_brace(self):
        formatter = SnakemakeFormatter()
        with pytest.raises(ValueError, match="unexpected '{' in field name"):
            list(formatter.parse("prefix_{subject{inner}suffix}"))

    def test_parse_doesnt_support_braces_in_constraints(self):
        formatter = SnakemakeFormatter()
        with pytest.raises(ValueError, match="unexpected '{' in field name"):
            list(formatter.parse("prefix_{subject,constraint{inner}suffix}"))

    def test_parse_raises_error_for_lone_closing_brace(self):
        formatter = SnakemakeFormatter()
        with pytest.raises(ValueError, match="unexpected '}' in string"):
            list(formatter.parse("prefix_}subject,constraint{inner}suffix}"))


class TestFormatting:
    """Basic tests for SnakemakeFormatter.format() method."""

    @given(
        literals=st.lists(literals()),
        field_names=st.lists(
            # Formatter splits field names on . and [ for indexing.
            field_names(exclude_characters=".[", min_size=1).filter(
                # Digits and empty strings are used to query args instead of kwargs
                # Values ending with _ or / are rejected by _validate
                lambda s: not s.isdigit() and not s.endswith(("/", "_"))
            )
        ),
    )
    def test_get_value_returns_kwargs_directly_when_present(
        self, literals: list[str], field_names: list[str]
    ):
        """Test that get_value() returns kwargs value directly when key exists."""
        wildcards = [f"{{{name},constraint}}" for name in field_names]
        path = "".join(itx.interleave_longest(literals, wildcards))
        target = "".join(
            itx.interleave_longest([s.format() for s in literals], field_names)
        )
        assert (
            SnakemakeFormatter().vformat(path, (), {n: n for n in field_names})
        ) == target

    @given(
        literals=st.lists(literals()),
        field_names=st.lists(st.integers(min_value=0, max_value=10**18)),
    )
    def test_numbered_format_fields_fallback_on_native_formatter(
        self, literals: list[str], field_names: list[int]
    ):
        """Test that get_value() returns kwargs value directly when key exists."""
        wildcards = [f"{{{name},constraint}}" for name in field_names]
        path = "".join(itx.interleave_longest(literals, wildcards))
        target = "".join(
            itx.interleave_longest(
                [s.format() for s in literals], map(str, field_names)
            )
        )
        assert (
            SnakemakeFormatter().vformat(path, range(max([*field_names, 0]) + 1), {})
        ) == target

    @given(
        literals=st.lists(literals()),
        field_names=st.lists(st.just("")),
    )
    def test_blank_format_fields_fallback_on_native_formatter(
        self, literals: list[str], field_names: list[str]
    ):
        """Test that get_value() returns kwargs value directly when key exists."""
        wildcards = [f"{{{name},constraint}}" for name in field_names]
        path = "".join(itx.interleave_longest(literals, wildcards))
        target = "".join(
            itx.interleave_longest(
                [s.format() for s in literals],
                map(str, range(len(field_names))),
            )
        )
        assert (
            SnakemakeFormatter().vformat(path, range(len(field_names) + 1), {})
        ) == target

    @given(
        literals=st.lists(literals()),
        field_names=st.lists(
            # Digits and empty strings are used to query args instead of kwargs
            safe_field_names(exclude_characters="_/", min_size=1).filter(
                lambda s: not s.isdigit()
            ),
            min_size=2,
        ),
    )
    def test_get_value_raises_error_for_missing_ordinary_entity(
        self, literals: list[str], field_names: list[str]
    ):
        unique = "".join(field_names)
        wildcards = [f"{{{name},constraint}}" for name in field_names] + [
            f"{{{unique},constraint}}"
        ]
        path = "".join(itx.interleave_longest(literals, wildcards))

        formatter = SnakemakeFormatter()
        with pytest.raises(KeyError, check=lambda e: e.args[0] == unique):
            formatter.vformat(path, (), {f: f for f in field_names})


class TestDirectoryWildcard:
    """Tests for directory wildcard functionality (_TAG_d_)."""

    @given(
        key=safe_field_names().map(lambda s: f"_{s}_d_"),
        value=st.text().map(lambda s: s + "/"),
    )
    def test_directory_returns_itself_when_given(self, key: str, value: str):
        formatter = SnakemakeFormatter()
        assert formatter.get_value(key, (), {key: value}) == value

    @given(entity=safe_field_names(), value=st.text())
    def test_directory_returns_default_when_corresponding_entity_given(
        self, entity: str, value: str
    ):
        assume(entity not in {"subject", "session"})
        key = f"_{entity}_d_"
        formatter = SnakemakeFormatter()
        target = f"{entity}-{value}/" if value else ""
        assert formatter.get_value(key, (), {entity: value}) == target

    @given(entity=safe_field_names(), value=st.text())
    def test_directory_raises_error_when_entity_missing(self, entity: str, value: str):
        key = f"_{entity}_d_"
        formatter = SnakemakeFormatter()
        with pytest.raises(
            KeyError,
            check=lambda e: e.args[0]
            == f"Missing required entity '{entity}' for wildcard '{key}'",
        ):
            formatter.get_value(key, (), {entity + "X": value})

    @pytest.mark.parametrize(
        ("template", "kwargs", "expected"),
        [
            ("{_subject_d_,constraint}anat", {"subject": "001"}, "sub-001/anat"),
            ("{_session_d_,constraint}func", {"session": "01"}, "ses-01/func"),
        ],
    )
    def test_directory_subject_and_session_use_sub_ses_tags(
        self, template: str, kwargs: dict[str, str], expected: str
    ):
        """Test directory wildcard uses 'sub' and 'ses' tags for subject/session."""
        formatter = SnakemakeFormatter()
        result = formatter.format(template, **kwargs)
        assert result == expected


class TestDummyWildcard:
    """Tests for dummy wildcard functionality (_TAG_)."""

    @given(
        key=safe_field_names().map(lambda s: f"_{s}_"),
        value=st.text().filter(lambda s: not s.endswith(("/", "_"))),
    )
    def test_dummy_returns_itself_when_given(self, key: str, value: str):
        """Test dummy wildcard returns assigned value when directly given as key."""
        formatter = SnakemakeFormatter()
        assert formatter.get_value(key, (), {key: value}) == value

    @given(
        entity=safe_field_names(),
        value=st.text(),
        underscore=st.sampled_from(["_", ""]),
    )
    def test_dummy_returns_default_when_corresponding_entity_given(
        self, entity: str, value: str, underscore: str
    ):
        assume(entity not in {"subject", "session"})
        key = f"_{entity}_"
        formatter = SnakemakeFormatter()
        formatter._underscore = underscore
        target = f"{underscore}{entity}-" if value else ""
        assert formatter.get_value(key, (), {entity: value}) == target

    @given(entity=safe_field_names(), value=st.text())
    def test_dummy_raises_error_when_entity_missing(self, entity: str, value: str):
        """Test dummy wildcard raises error when entity is missing."""
        key = f"_{entity}_"
        formatter = SnakemakeFormatter()
        with pytest.raises(
            KeyError,
            check=lambda e: e.args[0]
            == f"Missing required entity '{entity}' for wildcard '{key}'",
        ):
            formatter.get_value(key, (), {entity + "X": value})

    @pytest.mark.parametrize(
        ("template", "kwargs", "expected"),
        [
            ("{_subject_,constraint}", {"subject": "001"}, "sub-"),
            ("{_session_,constraint}", {"session": "01"}, "ses-"),
        ],
    )
    def test_dummy_subject_and_session_use_sub_ses_tags(
        self, template: str, kwargs: dict[str, str], expected: str
    ):
        """Test dummy wildcard uses 'sub' and 'ses' tags for subject/session."""
        formatter = SnakemakeFormatter()
        result = formatter.format(template, **kwargs)
        assert result == expected

    def test_single_underscore_treated_as_ordinary_wildcard(self):
        formatter = SnakemakeFormatter()
        assert formatter.format("{_}", _="foo") == "foo"
        with pytest.raises(KeyError, match="^'_'$"):
            formatter.format("{_}")


class TestUnderscoreSquelching:
    """Test all formats twice to ensure vformat properly resets state."""

    def test_formatter_state_resets_each_vformat(self):
        formatter = SnakemakeFormatter()
        assert formatter.vformat("{_foo_}", (), {"foo": "foo"}) == "foo-"
        assert formatter.vformat("{_foo_}", (), {"foo": "foo"}) == "foo-"

    def test_formatter_state_resets_each_format(self):
        formatter = SnakemakeFormatter()
        assert formatter.format("{_foo_}", foo="foo") == "foo-"
        assert formatter.format("{_foo_}", foo="foo") == "foo-"

    @example(literal="/", entity="foo")
    @example(literal="_", entity="foo")
    @given(literal=literals(), entity=safe_field_names() | st.just("_"))
    def test_squelching_occurs_at_beginning_of_string_and_after_slash(
        self, literal: str, entity: str
    ):
        assume(entity not in {"subject", "session"})
        key = f"_{entity}_"
        formatter = SnakemakeFormatter()
        target = literal.format() + (
            "_"
            if literal and literal[-1] not in SnakemakeFormatter.UNDERSCORE_SQUELCHERS
            else ""
        )
        if entity != "_":
            target += entity + "-"
        assert formatter.vformat(f"{literal}{{{key}}}", (), {entity: "foo"}) == target

    @given(
        literal=literals().map(lambda s: s.rstrip("/_")),
        value=st.text().filter(lambda s: not s.endswith(("/", "_"))),
        entity=st.sampled_from(["_", "entity"]),
        missing=st.booleans(),
    )
    def test_squelching_after_ordinary_wildcard(
        self, literal: str, value: str, entity: str, missing: bool
    ):
        formatter = SnakemakeFormatter(allow_missing=missing)

        if missing:
            underscore = "_" if literal else str(SnakemakeWildcards.underscore)
        elif value:
            underscore = "_"
        else:
            underscore = ""

        target = (
            literal.format()
            + ("_" if literal and not value and not missing else "")
            + ("{before}" if missing else value)
            + underscore
            + (entity + "-" if entity != "_" else "")
        )

        kwargs = {entity: "foo"}
        if not missing:
            kwargs["before"] = value
        assert (
            formatter.format(f"{literal}{{before}}{{_{entity}_}}", **kwargs) == target
        )

    @given(
        literal=literals().map(lambda s: s.rstrip("/_")),
        before=st.sampled_from(["_", "before"]),
        skip=st.booleans(),
        entity=st.sampled_from(["_", "entity"]),
    )
    def test_squelching_after_dummy_wildcard(
        self, literal: str, before: str, skip: bool, entity: str
    ):
        formatter = SnakemakeFormatter()
        target = (
            literal.format()
            + ("_" if literal else "")
            + (f"{before}-_" if not skip and before != "_" else "")
            + (entity + "-" if entity != "_" else "")
        )

        kwargs = {entity: "foo", before: "" if skip else "foo"}
        assert (
            formatter.format(f"{literal}{{_{before}_}}{{_{entity}_}}", **kwargs)
            == target
        )

    @given(
        literal=literals().map(lambda s: s.rstrip("/_")),
        entity=st.sampled_from(["_", "entity"]),
    )
    def test_squelching_after_missing_dummy_wildcard(self, literal: str, entity: str):
        formatter = SnakemakeFormatter(allow_missing=True)
        formatted_before = "{_before_}" + (
            str(SnakemakeWildcards.underscore) if not literal else "_"
        )
        target = (
            literal.format()
            + formatted_before
            + (entity + "-" if entity != "_" else "")
        )

        kwargs = {entity: "foo"}
        assert (
            formatter.format(f"{literal}{{_before_}}{{_{entity}_}}", **kwargs) == target
        )

    @pytest.mark.parametrize(
        ("before_wildcard", "before_entity", "before_value"),
        [
            ("__d__", "datatype", "/"),
            ("_dir_d_", "dir", "dir-foo/"),
        ],
    )
    @given(
        literal=literals(),
        skip=st.booleans(),
        entity=st.sampled_from(["_", "entity"]),
    )
    def test_squelching_occurs_after_slash_wildcard(
        self,
        literal: str,
        skip: bool,
        entity: str,
        before_wildcard: str,
        before_entity: str,
        before_value: str,
    ):
        formatter = SnakemakeFormatter()
        target = (
            literal.format()
            + (
                before_value
                if not skip
                else "_"
                if literal and not literal.endswith(("/", "_"))
                else ""
            )
            + (entity + "-" if entity != "_" else "")
        )
        assert (
            formatter.vformat(
                f"{literal}{{{before_wildcard}}}{{_{entity}_}}",
                (),
                {before_entity: "" if skip else "foo", entity: "foo"},
            )
            == target
        )

    @given(
        literal=literals(),
        entity=st.sampled_from(["_", "entity"]),
        before_wildcard=st.sampled_from(["__d__", "_dir_d_"]),
    )
    def test_squelching_occurs_after_missing_slash_wildcard(
        self,
        literal: str,
        entity: str,
        before_wildcard: str,
    ):
        formatter = SnakemakeFormatter(allow_missing=True)
        target = (
            literal.format()
            + f"{{{before_wildcard}}}"
            + (entity + "-" if entity != "_" else "")
        )
        assert (
            formatter.format(
                f"{literal}{{{before_wildcard}}}{{_{entity}_}}", **{entity: "foo"}
            )
            == target
        )

    @given(
        before=st.sampled_from(
            [
                "__d__",
                "datatype",
                "wildcard",
                "_wildcard_",
                "_wildcard_d_",
                "_empty_",
                "_empty_d_",
                "_empty_d_",
            ]
        ).map(lambda s: f"{{{s}}}")
        | literals().map(lambda s: s.rstrip("_")),
        wildcard=st.sampled_from(["__d__", "_dir_d_"]),
    )
    def test_never_underscore_before_slash_wildcard(self, before: str, wildcard: str):
        formatter = SnakemakeFormatter()
        result = formatter.vformat(
            f"{before}{{{wildcard}}}",
            (),
            {"datatype": "foo", "wildcard": "foo", "empty": "", "dir": "foo"},
        )
        tail = "/" if wildcard == "__d__" else "dir-foo/"
        assert result.endswith(tail)
        assert not result.endswith(f"_{tail}")


class TestDWildcard:
    """Tests for __d__ special wildcard."""

    @pytest.mark.parametrize(("datatype", "expected"), [("anat", "/"), ("", "")])
    def test_d_returns_slash_when_datatype_present(self, datatype: str, expected: str):
        """Test __d__ returns / when datatype is present and non-blank."""
        formatter = SnakemakeFormatter()
        assert formatter.get_value("__d__", (), {"datatype": datatype}) == expected

    def test_d_raises_error_when_datatype_missing(self):
        """Test __d__ raises error when datatype is missing."""
        formatter = SnakemakeFormatter()
        with pytest.raises(
            KeyError, match="Missing required entity 'datatype' for wildcard '__d__'"
        ):
            formatter.get_value("__d__", (), {})


class TestTripleUnderscore:
    """Tests for ___ special wildcard."""

    @pytest.mark.parametrize("underscore", ["_", ""])
    def test_triple_underscore_returns_based_on_squelch_state(self, underscore: str):
        formatter = SnakemakeFormatter()
        formatter._underscore = underscore
        result = formatter.get_value("___", (), {})
        assert result == underscore


class TestIntegration:
    """Integration tests for complete formatting scenarios."""

    def test_format_with_required_entities_only(self):
        """Test formatting with only required entities."""
        formatter = SnakemakeFormatter()
        result = formatter.format(
            "sub-{subject}/anat/sub-{subject}_T1w.nii.gz", subject="001"
        )
        assert result == "sub-001/anat/sub-001_T1w.nii.gz"

    def test_format_with_optional_entities(self):
        """Test formatting with optional entities."""
        formatter = SnakemakeFormatter()
        result = formatter.format(
            "sub-{subject}{_acq_}{acq}{_run_}{run}_T1w.nii.gz",
            subject="001",
            acq="mprage",
            run="01",
        )
        assert result == "sub-001_acq-mprage_run-01_T1w.nii.gz"

    def test_format_with_directory_wildcards(self):
        """Test formatting with directory wildcards."""
        formatter = SnakemakeFormatter()
        result = formatter.format(
            "{_subject_d_}{_session_d_}anat/sub-{subject}{_session_}{session}"
            "_T1w.nii.gz",
            subject="001",
            session="01",
        )
        assert result == "sub-001/ses-01/anat/sub-001_ses-01_T1w.nii.gz"

    def test_format_with_subject_session_mapping(self):
        """Test subject and session mapping to sub and ses."""
        formatter = SnakemakeFormatter()
        result = formatter.format(
            "{_subject_d_}{_session_d_}func/{_subject_}{subject}{_session_}{session}"
            "{___}task-rest_bold.nii.gz",
            subject="002",
            session="02",
        )
        assert result == "sub-002/ses-02/func/sub-002_ses-02_task-rest_bold.nii.gz"

    def test_format_validates_against_discussion_examples(self):
        """Test formatting against various real-world BIDS patterns."""
        formatter = SnakemakeFormatter()

        # Example 1: Simple BIDS path with subject and datatype
        result = formatter.format(
            "bids/{_subject_d_}{datatype}{__d__}{_subject_}{subject}_T1w{extension}",
            datatype="anat",
            subject="001",
            extension=".nii.gz",
        )
        assert result == "bids/sub-001/anat/sub-001_T1w.nii.gz"

        # Example 2: Optional acquisition and run
        result = formatter.format(
            "sub-{subject}{_acq_}{acq}{_run_}{run}_bold{extension}",
            subject="001",
            acq="",
            run="01",
            extension=".nii.gz",
        )
        assert result == "sub-001_run-01_bold.nii.gz"

        # Example 3: Full path with all components
        result = formatter.format(
            "{_subject_d_}{_session_d_}{datatype}{__d__}{_subject_}{subject}{_session_}"
            "{session}{_acq_}{acq}{___}{suffix}{extension}",
            subject="003",
            session="01",
            datatype="func",
            acq="bold",
            suffix="rest",
            extension=".nii.gz",
        )
        assert result == "sub-003/ses-01/func/sub-003_ses-01_acq-bold_rest.nii.gz"

    def test_format_with_optional_entities_missing(self):
        """Test that optional entities can be omitted (blank or None)."""
        formatter = SnakemakeFormatter()

        # With blank strings
        result = formatter.format(
            "{_subject_}{subject}{_acq_}{acq}{_run_}{run}{___}bold.nii.gz",
            subject="001",
            acq="",
            run="",
        )
        assert result == "sub-001_bold.nii.gz"

    def test_format_complex_path_with_constraints(self):
        """Test formatting with wildcard constraints (should be stripped)."""
        formatter = SnakemakeFormatter()
        result = formatter.format(
            "sub-{subject,\\d+}{_acq_,\\w+}{acq}_T1w.nii.gz",
            subject="001",
            acq="mprage",
        )
        assert result == "sub-001_acq-mprage_T1w.nii.gz"


class TestValidation:
    """Tests for _validate error conditions in SnakemakeFormatter."""

    @given(value=st.text(min_size=1).filter(lambda s: s != "/"))
    def test_d_wildcard_rejects_non_slash_value(self, value: str):
        """__d__ assigned a value other than '/' raises ValueError."""
        formatter = SnakemakeFormatter()
        with pytest.raises(ValueError, match="__d__ assigned invalid value"):
            formatter.format("{__d__}", __d__=value)

    @pytest.mark.parametrize("value", ["/", ""])
    def test_d_wildcard_accepts_valid_values(self, value: str):
        """__d__ assigned '/' or '' passes validation."""
        formatter = SnakemakeFormatter()
        assert formatter.format("{__d__}", __d__=value) == value

    @given(value=st.text(min_size=1).filter(lambda s: s != "_"))
    def test_triple_underscore_rejects_non_underscore_value(self, value: str):
        """___ assigned a value other than '_' raises ValueError."""
        formatter = SnakemakeFormatter()
        with pytest.raises(ValueError, match="___ assigned invalid value"):
            formatter.format("{___}", ___=value)

    @pytest.mark.parametrize("value", ["_", ""])
    def test_triple_underscore_accepts_valid_values(self, value: str):
        """___ assigned '_' or '' passes validation."""
        formatter = SnakemakeFormatter()
        assert formatter.format("{___}", ___=value) == value

    @given(value=st.text(min_size=1).filter(lambda s: not s.endswith("/")))
    def test_directory_wildcard_rejects_value_not_ending_in_slash(self, value: str):
        """Directory wildcard assigned value not ending in '/' raises ValueError."""
        formatter = SnakemakeFormatter()
        with pytest.raises(ValueError, match="does not end in '/'"):
            formatter.get_value("_entity_d_", (), {"_entity_d_": value})

    @given(value=st.text().map(lambda s: s + "/"))
    def test_directory_wildcard_accepts_value_ending_in_slash(self, value: str):
        """Directory wildcard assigned value ending in '/' passes validation."""
        formatter = SnakemakeFormatter()
        assert formatter.get_value("_entity_d_", (), {"_entity_d_": value}) == value

    @given(
        bad_suffix=st.sampled_from(["_", "/"]),
        value=st.text(),
    )
    def test_ordinary_entity_rejects_value_ending_with_disallowed_char(
        self, bad_suffix: str, value: str
    ):
        """Ordinary entity value ending with '_' or '/' raises ValueError."""
        formatter = SnakemakeFormatter()
        with pytest.raises(ValueError, match="ending with disallowed character"):
            formatter.format("{entity}", entity=value + bad_suffix)

    @given(
        value=st.text(min_size=1).filter(lambda s: not s.endswith(("/", "_"))),
    )
    def test_ordinary_entity_accepts_valid_value(self, value: str):
        """Ordinary entity value not ending with '_' or '/' passes validation."""
        formatter = SnakemakeFormatter()
        assert formatter.format("{entity}", entity=value) == value

    @given(
        bad_suffix=st.sampled_from(["_", "/"]),
        value=st.text(),
    )
    def test_integer_key_validates_value(self, bad_suffix: str, value: str):
        """Integer key with value ending in disallowed char raises ValueError."""
        formatter = SnakemakeFormatter()
        with pytest.raises(ValueError, match="ending with disallowed character"):
            formatter.vformat("{0}", [value + bad_suffix], {})

    @given(
        value=st.text(min_size=1).filter(lambda s: not s.endswith(("/", "_"))),
    )
    def test_integer_key_accepts_valid_value(self, value: str):
        """Integer key with normal value passes validation."""
        formatter = SnakemakeFormatter()
        assert formatter.vformat("{0}", [value], {}) == value


class TestAllowMissing:
    """Tests for allow_missing behavior in SnakemakeFormatter."""

    # MARK
    @example(name="run", constraint=r"\d+")
    @example(name="__d__", constraint=r"\w+")
    @example(name="_dir_d_", constraint="")
    @example(name="_acq_", constraint="")
    @given(
        name=safe_field_names(min_size=1).filter(lambda s: not s.isdigit()),
        constraint=constraints() | st.just(""),
    )
    def test_allow_missing_preserves_wildcard_from_template(
        self, name: str, constraint: str
    ):
        """Missing wildcard (with or without constraint) is preserved verbatim."""
        formatter = SnakemakeFormatter(allow_missing=True)
        template = f"{{{name},{constraint}}}" if constraint else f"{{{name}}}"
        result = formatter.format(template)
        assert result == template

    def test_allow_missing_resets_between_vformat_calls(self):
        """_field_constraints is reset on each vformat call."""
        formatter = SnakemakeFormatter(allow_missing=True)
        # First call with constraint
        r1 = formatter.format(r"sub-{subject,\d+}")
        assert r1 == r"sub-{subject,\d+}"
        # Second call without constraint: should use the new (no-constraint) form
        r2 = formatter.format("{subject}")
        assert r2 == "{subject}"

    @given(
        name=safe_field_names(min_size=1).filter(lambda s: not s.isdigit()),
        constraint_a=constraints() | st.just(""),
        constraint_b=constraints() | st.just(""),
    )
    def test_allow_missing_preserves_each_fields_constraint(
        self, name: str, constraint_a: str, constraint_b: str
    ):
        """Same wildcard with different constraints preserves each independently."""
        formatter = SnakemakeFormatter(allow_missing=True)
        field_a = f"{{{name},{constraint_a}}}" if constraint_a else f"{{{name}}}"
        field_b = f"{{{name},{constraint_b}}}" if constraint_b else f"{{{name}}}"
        template = f"{field_a}_lit_{field_b}"
        result = formatter.format(template)
        assert result == template


@pytest.mark.skipif(not _HAS_RUST_PARSE, reason="Rust extension not built")
class TestRustParityParse:
    """Verify that the Rust-backed parse() is byte-for-byte identical to the
    pure-Python implementation across a range of inputs.

    Each helper forces a specific code path so regressions are easy to spot.

    Note: every test here calls ``_parse_python()`` directly, so the
    pure-Python fallback is also exercised.  When the Rust extension is absent
    (e.g. in a plain ``pip install`` environment) the existing ``TestParse``
    class covers ``parse()`` via the Python path.
    """

    @staticmethod
    def _both(template: str) -> tuple[list, list]:
        """Return (python_results, rust_results) for *template*."""
        py = SnakemakeFormatter()
        rust = SnakemakeFormatter()

        py_entries = list(py._parse_python(template))
        rust_entries = list(rust._parse_rust(template))
        return py_entries, rust_entries

    # ---- yield-tuple parity ----------------------------------------------

    def test_parity_pure_literal(self):
        py, rs = self._both("just a literal")
        assert py == rs

    def test_parity_simple_field(self):
        py, rs = self._both("{field}")
        assert py == rs

    def test_parity_literal_then_field(self):
        py, rs = self._both("prefix_{field}")
        assert py == rs

    def test_parity_field_then_literal(self):
        py, rs = self._both("{field}_suffix")
        assert py == rs

    def test_parity_doubled_open_brace(self):
        py, rs = self._both("a{{b")
        assert py == rs

    def test_parity_doubled_close_brace(self):
        py, rs = self._both("a}}b")
        assert py == rs

    def test_parity_constraint_field(self):
        py, rs = self._both(r"{subject,\d+}")
        assert py == rs

    def test_parity_constraint_field_with_literal(self):
        py, rs = self._both(r"sub-{subject,\d+}_T1w")
        assert py == rs

    def test_parity_multiple_fields(self):
        py, rs = self._both("{a}_{b}_{c}")
        assert py == rs

    def test_parity_mixed_constraints_and_plain(self):
        py, rs = self._both(r"{a,\w+}_{b}_{c,\d+}")
        assert py == rs

    def test_parity_empty_string(self):
        py, rs = self._both("")
        assert py == rs

    # ---- error parity ----------------------------------------------------

    @pytest.mark.parametrize(
        "bad",
        [
            "prefix_{subject",
            "prefix_{",
            "prefix_}",
            "prefix_{subject{inner}suffix}",
            "prefix_{subject,constraint{inner}suffix}",
        ],
    )
    def test_both_raise_same_error(self, bad: str):
        with pytest.raises(ValueError) as py_exc:
            list(SnakemakeFormatter()._parse_python(bad))
        with pytest.raises(ValueError) as rs_exc:
            list(SnakemakeFormatter()._parse_rust(bad))
        assert py_exc.value.args == rs_exc.value.args

    # ---- side-effect parity: _underscore ---------------------------------

    def test_underscore_after_pure_literal(self):
        py = SnakemakeFormatter()
        rs = SnakemakeFormatter()
        list(py._parse_python("abc"))
        list(rs._parse_rust("abc"))
        assert py._underscore == rs._underscore

    def test_underscore_after_doubled_brace(self):
        py = SnakemakeFormatter()
        rs = SnakemakeFormatter()
        list(py._parse_python("{{"))
        list(rs._parse_rust("{{"))
        assert py._underscore == rs._underscore

    def test_underscore_after_slash_literal(self):
        py = SnakemakeFormatter()
        rs = SnakemakeFormatter()
        list(py._parse_python("a/_{field}"))
        list(rs._parse_rust("a/_{field}"))
        assert py._underscore == rs._underscore

    def test_underscore_after_underscore_literal(self):
        py = SnakemakeFormatter()
        rs = SnakemakeFormatter()
        list(py._parse_python("a__{field}"))
        list(rs._parse_rust("a__{field}"))
        assert py._underscore == rs._underscore

    # ---- side-effect parity: _current_field ------------------------------

    def test_current_field_with_constraint(self):
        py = SnakemakeFormatter()
        rs = SnakemakeFormatter()
        list(py._parse_python(r"{field,\d+}"))
        list(rs._parse_rust(r"{field,\d+}"))
        assert py._current_field == rs._current_field

    def test_current_field_without_constraint(self):
        # With no constraint, _parse_rust sets _current_field to the field name
        # (field_name + "" == field_name), whereas _parse_python sets it to None.
        # Both produce identical formatted output since `_current_field or key` gives
        # the same result in both cases.
        rs = SnakemakeFormatter()
        list(rs._parse_rust("{field}"))
        assert rs._current_field == "field"

    # ---- allow_missing parity (exercises _current_field in get_value) ---

    @example(name="run", constraint=r"\d+")
    @given(
        name=safe_field_names(min_size=1).filter(lambda s: not s.isdigit()),
        constraint=constraints() | st.just(""),
    )
    def test_allow_missing_identical_output(self, name: str, constraint: str):
        template = f"{{{name},{constraint}}}" if constraint else f"{{{name}}}"
        py = SnakemakeFormatter(allow_missing=True)
        rs = SnakemakeFormatter(allow_missing=True)
        # Override to force each path regardless of _HAS_RUST_PARSE
        py_result = list(py._parse_python(template))
        rs_result = list(rs._parse_rust(template))
        assert py_result == rs_result

    # ---- broad property-based parity -------------------------------------

    @example(literals=["}}"], wildcards=["{}", "{}"])
    @given(
        literals=st.lists(literals()),
        wildcards=st.lists(
            field_names(exclude_characters="!:").map(lambda s: f"{{{s}}}")
        ),
    )
    def test_rust_parse_matches_python_parse(
        self, literals: list[str], wildcards: list[str]
    ):
        path = "".join(itx.interleave_longest(literals, wildcards))
        py, rs = self._both(path)
        assert py == rs


def _bids_args():
    return st.dictionaries(
        keys=sb_st.bids_entity().map(lambda e: e.wildcard),
        values=st.text(st.characters(exclude_characters="_-./{}"), min_size=1),
    ).filter(lambda s: set(s) - {"datatype", "extension"})


@given(bidsargs=_bids_args())
def test_single_and_multistep_format_equivalent(bidsargs: dict[str, str]):
    template = bids("", **dict.fromkeys(bidsargs, OPTIONAL_WILDCARD))
    formatter = SnakemakeFormatter(allow_missing=True)
    single_step = formatter.format(template, **bidsargs)
    multistep = template
    for k, v in bidsargs.items():
        multistep = formatter.format(multistep, **{k: v})
    assert single_step == multistep
