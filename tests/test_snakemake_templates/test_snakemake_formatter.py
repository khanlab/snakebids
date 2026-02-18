"""Tests for snakemake template formatting utilities."""

from __future__ import annotations

import string

import hypothesis.strategies as st
import more_itertools as itx
import pytest
from hypothesis import assume, example, given

from snakebids.utils.snakemake_templates import SnakemakeFormatter
from tests.helpers import Benchmark
from tests.test_snakemake_templates.strategies import (
    constraints,
    field_names,
    literals,
    safe_field_names,
)


class TestParserBenchmarks:
    times = 10000
    text = "fewaigq{{eafwefao{wafwa,waffwahaw,a:fa}{fea!f:wagwa}jawjigaw{fwa,afewaew}k"

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


class TestParse:
    """Tests for SnakemakeFormatter.parse() method."""

    @given(
        literals=st.lists(literals().map(lambda s: s.replace("}", ""))),
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
        """Test that raises ValueError for a trailing opening brace."""
        formatter = SnakemakeFormatter()
        with pytest.raises(ValueError, match="expected '}' before end of string"):
            list(formatter.parse("prefix_{"))

    def test_parse_raises_error_for_nested_opening_brace(self):
        """Test that raises ValueError for nested opening brace in field name."""
        formatter = SnakemakeFormatter()
        with pytest.raises(ValueError, match="unexpected '{' in field name"):
            list(formatter.parse("prefix_{subject{inner}suffix}"))

    def test_parse_doesnt_support_braces_in_constraints(self):
        formatter = SnakemakeFormatter()
        with pytest.raises(ValueError, match="unexpected '{' in field name"):
            list(formatter.parse("prefix_{subject,constraint{inner}suffix}"))


class TestFormatting:
    """Basic tests for SnakemakeFormatter.format() method."""

    @given(
        literals=st.lists(literals()),
        field_names=st.lists(
            # Formatter splits field names on . and [ for indexing.
            field_names(exclude_characters=".[", min_size=1).filter(
                # Digits and empty strings are used to query args instead of kwargs
                lambda s: not s.isdigit()
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
            itx.interleave_longest(
                [s.replace("{{", "{") for s in literals], field_names
            )
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
                [s.replace("{{", "{") for s in literals], map(str, field_names)
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
                [s.replace("{{", "{") for s in literals],
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
            safe_field_names(min_size=1).filter(lambda s: not s.isdigit()),
            min_size=2,
        ),
    )
    def test_get_value_raises_error_for_missing_ordinary_entity(
        self, literals: list[str], field_names: list[str]
    ):
        """Test that get_value() returns kwargs value directly when key exists."""
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

    @given(key=safe_field_names().map(lambda s: f"_{s}_d_"), value=st.text())
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

    @given(key=safe_field_names().map(lambda s: f"_{s}_"), value=st.text())
    def test_dummy_returns_itself_when_given(self, key: str, value: str):
        """Test dummy wildcard returns assigned value when directly given as key."""
        formatter = SnakemakeFormatter()
        assert formatter.get_value(key, (), {key: value}) == value

    @pytest.mark.parametrize("squelch_underscore", [True, False])
    @given(entity=safe_field_names(), value=st.text())
    def test_dummy_returns_default_when_corresponding_entity_given(
        self, entity: str, value: str, squelch_underscore: bool
    ):
        """Test dummy wildcard returns the correct format based on squelch state."""
        assume(entity not in {"subject", "session"})
        key = f"_{entity}_"
        formatter = SnakemakeFormatter()
        formatter.squelch_underscore = squelch_underscore
        if value:
            target = f"{entity}-" if squelch_underscore else f"_{entity}-"
        else:
            target = ""
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
        target = literal.replace("{{", "{") + (
            "_"
            if literal and literal[-1] not in SnakemakeFormatter.UNDERSCORE_SQUELCHERS
            else ""
        )
        if entity != "_":
            target += entity + "-"
        assert formatter.vformat(f"{literal}{{{key}}}", (), {entity: "foo"}) == target

    @given(
        literal=literals().map(lambda s: s.rstrip("/_")),
        value=st.text() | st.sampled_from(["/", "_"]),
        entity=st.sampled_from(["_", "entity"]),
    )
    def test_squelching_after_ordinary_wildcard(
        self, literal: str, value: str, entity: str
    ):
        formatter = SnakemakeFormatter()
        target = f"{entity}-" if entity != "_" else ""
        if value and value[-1] not in SnakemakeFormatter.UNDERSCORE_SQUELCHERS:
            target = f"{value}_{target}"
        elif literal and not value:
            target = f"_{target}"
        else:
            target = f"{value}{target}"

        target = literal.replace("{{", "{") + f"{target}"

        assert (
            formatter.vformat(
                f"{literal}{{before}}{{_{entity}_}}",
                (),
                {"before": value, entity: "foo"},
            )
            == target
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
            literal.replace("{{", "{")
            + ("_" if literal else "")
            + (f"{before}-_" if not skip and before != "_" else "")
            + (entity + "-" if entity != "_" else "")
        )

        assert (
            formatter.vformat(
                f"{literal}{{_{before}_}}{{_{entity}_}}",
                (),
                {before: "" if skip else "foo", entity: "foo"},
            )
            == target
        )

    @pytest.mark.parametrize(
        ("before_wildcard", "before_entity", "before_value"),
        [
            ("__d__", "datatype", "/"),
            ("_dir_d_", "dir", "dir-foo/"),
        ],
    )
    @given(
        literal=literals().map(lambda s: s.rstrip("/_")),
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
            literal.replace("{{", "{")
            + (before_value if not skip else "_" if literal else "")
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

    @pytest.mark.parametrize(
        ("squelch_underscore", "expected"),
        [
            (False, "_"),
            (True, ""),
        ],
    )
    def test_triple_underscore_returns_based_on_squelch_state(
        self, squelch_underscore: bool, expected: str
    ):
        """Test ___ returns underscore when not squelched, blank when squelched."""
        formatter = SnakemakeFormatter()
        formatter.squelch_underscore = squelch_underscore
        result = formatter.get_value("___", (), {})
        assert result == expected


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
