from __future__ import annotations

import logging
import shutil
from pathlib import Path

import pandas as pd
import pytest
from upath import UPath

from snakebids.core.input_generation import generate_inputs
from snakebids.core.participants import query_participants, read_participants_tsv
from snakebids.exceptions import ConfigError
from snakebids.types import InputsConfig

PARTICIPANTS = (
    "participant_id\tage\tgroup\n"
    "sub-001\t25\tcontrol\n"
    "sub-002\t67\tpatient\n"
    "sub-003\tn/a\tpatient\n"
)

PYBIDS_INPUTS: InputsConfig = {
    "t1": {"filters": {"suffix": "T1w"}, "wildcards": ["subject", "acquisition"]}
}


def write_participants(root: Path | UPath, contents: str = PARTICIPANTS) -> None:
    (root / "participants.tsv").write_text(contents)


@pytest.fixture
def bids_dir(tmp_path: Path) -> Path:
    root = tmp_path / "bids"
    shutil.copytree("tests/data/bids_t1w", root)
    (root / "dataset_description.json").write_text(
        '{"Name": "test", "BIDSVersion": "1.8.0"}'
    )
    return root


class TestReadParticipantsTsv:
    def test_returns_none_when_missing(self, tmp_path: Path):
        assert read_participants_tsv(tmp_path) is None

    def test_indexed_by_subject_label(self, tmp_path: Path):
        write_participants(tmp_path)
        table = read_participants_tsv(tmp_path)
        assert table is not None
        assert table.index.name == "subject"
        assert list(table.index) == ["001", "002", "003"]
        assert list(table["participant_id"]) == ["sub-001", "sub-002", "sub-003"]

    def test_na_parsed_and_numeric_columns_inferred(self, tmp_path: Path):
        write_participants(tmp_path)
        table = read_participants_tsv(tmp_path)
        assert table is not None
        assert pd.api.types.is_numeric_dtype(table["age"])
        assert pd.isna(table.loc["003", "age"])

    def test_numeric_looking_labels_kept_as_strings(self, tmp_path: Path):
        write_participants(tmp_path, "participant_id\tage\n001\t20\n")
        table = read_participants_tsv(tmp_path)
        assert table is not None
        assert list(table.index) == ["001"]

    def test_missing_participant_id_raises_error(self, tmp_path: Path):
        write_participants(tmp_path, "subject\tage\nsub-001\t20\n")
        with pytest.raises(ConfigError, match="participant_id"):
            read_participants_tsv(tmp_path)

    def test_unparseable_file_raises_error(self, tmp_path: Path):
        write_participants(tmp_path, "")
        with pytest.raises(ConfigError, match="Could not parse"):
            read_participants_tsv(tmp_path)

    def test_reads_remote_dataset(self):
        # fsspec's in-memory filesystem stands in for gs://, s3://, etc.
        root = UPath("memory://bucket/remote_ds")
        root.mkdir(parents=True, exist_ok=True)
        write_participants(root)
        table = read_participants_tsv("memory://bucket/remote_ds")
        assert table is not None
        assert list(table.index) == ["001", "002", "003"]

    def test_returns_none_when_missing_from_remote_dataset(self):
        assert read_participants_tsv("memory://bucket/empty_ds") is None

    def test_duplicate_ids_warn(self, tmp_path: Path, caplog: pytest.LogCaptureFixture):
        write_participants(tmp_path, "participant_id\nsub-001\nsub-001\n")
        with caplog.at_level(logging.WARNING):
            read_participants_tsv(tmp_path)
        assert "duplicate" in caplog.text


class TestQueryParticipants:
    @pytest.fixture
    def table(self, tmp_path: Path) -> pd.DataFrame:
        write_participants(tmp_path)
        table = read_participants_tsv(tmp_path)
        assert table is not None
        return table

    def test_numeric_query(self, table: pd.DataFrame):
        assert query_participants(table, "age > 30") == ["002"]

    def test_combined_query(self, table: pd.DataFrame):
        assert query_participants(table, "group == 'patient' and age < 30") == []
        assert query_participants(table, "group == 'patient'") == ["002", "003"]

    def test_query_on_subject_index(self, table: pd.DataFrame):
        assert query_participants(table, "subject in ['001', '003']") == [
            "001",
            "003",
        ]

    def test_unknown_column_raises_error_listing_columns(self, table: pd.DataFrame):
        with pytest.raises(ConfigError, match="Available columns: participant_id"):
            query_participants(table, "sex == 'F'")

    def test_invalid_syntax_raises_error(self, table: pd.DataFrame):
        with pytest.raises(ConfigError, match="Could not evaluate"):
            query_participants(table, "age >")


class TestGenerateInputs:
    def test_participants_attached_to_dataset(self, bids_dir: Path):
        write_participants(bids_dir)
        dataset = generate_inputs(bids_dir, PYBIDS_INPUTS)
        assert dataset.participants is not None
        assert list(dataset.participants.index) == ["001", "002", "003"]

    def test_malformed_tsv_warns_without_filter(
        self, bids_dir: Path, caplog: pytest.LogCaptureFixture
    ):
        write_participants(bids_dir, "subject\tage\nsub-001\t20\n")
        with caplog.at_level(logging.WARNING):
            dataset = generate_inputs(bids_dir, PYBIDS_INPUTS)
        assert dataset.participants is None
        assert "participant_id" in caplog.text
        assert sorted(dataset.subjects) == ["001", "002"]

    def test_malformed_tsv_raises_with_filter(self, bids_dir: Path):
        write_participants(bids_dir, "subject\tage\nsub-001\t20\n")
        with pytest.raises(ConfigError, match="participant_id"):
            generate_inputs(bids_dir, PYBIDS_INPUTS, participant_filter="age > 0")

    def test_participants_none_without_tsv(self, bids_dir: Path):
        dataset = generate_inputs(bids_dir, PYBIDS_INPUTS)
        assert dataset.participants is None

    def test_filter_selects_subjects(self, bids_dir: Path):
        write_participants(bids_dir)
        dataset = generate_inputs(
            bids_dir, PYBIDS_INPUTS, participant_filter="age > 30"
        )
        assert dataset.subjects == ["002"]
        # table is not filtered
        assert dataset.participants is not None
        assert list(dataset.participants.index) == ["001", "002", "003"]

    def test_filter_intersects_with_participant_label(self, bids_dir: Path):
        write_participants(bids_dir)
        dataset = generate_inputs(
            bids_dir,
            PYBIDS_INPUTS,
            participant_filter="age > 0",
            participant_label=["002"],
        )
        assert dataset.subjects == ["002"]
        dataset = generate_inputs(
            bids_dir,
            PYBIDS_INPUTS,
            participant_filter="age > 30",
            participant_label="001",
        )
        assert dataset.subjects == []

    def test_filter_combines_with_exclude_participant_label(self, bids_dir: Path):
        write_participants(bids_dir)
        dataset = generate_inputs(
            bids_dir,
            PYBIDS_INPUTS,
            participant_filter="age > 0",
            exclude_participant_label="002",
        )
        assert dataset.subjects == ["001"]

    def test_unlisted_subjects_excluded_with_warning(
        self, bids_dir: Path, caplog: pytest.LogCaptureFixture
    ):
        write_participants(bids_dir, "participant_id\tage\nsub-002\t50\n")
        with caplog.at_level(logging.WARNING):
            dataset = generate_inputs(
                bids_dir, PYBIDS_INPUTS, participant_filter="age > 0"
            )
        assert dataset.subjects == ["002"]
        assert "not listed in participants.tsv" in caplog.text
        assert "001" in caplog.text

    def test_no_match_warns(self, bids_dir: Path, caplog: pytest.LogCaptureFixture):
        write_participants(bids_dir)
        with caplog.at_level(logging.WARNING):
            dataset = generate_inputs(
                bids_dir, PYBIDS_INPUTS, participant_filter="age > 100"
            )
        assert dataset.subjects == []
        assert "No participants matched" in caplog.text

    def test_filter_without_tsv_raises_error(self, bids_dir: Path):
        with pytest.raises(ConfigError, match="has no participants.tsv"):
            generate_inputs(bids_dir, PYBIDS_INPUTS, participant_filter="age > 30")

    def test_filter_applies_to_custom_paths(self, bids_dir: Path):
        write_participants(bids_dir)
        dataset = generate_inputs(
            bids_dir,
            {
                "t1": {
                    "filters": {},
                    "wildcards": ["subject"],
                    "custom_path": str(
                        bids_dir
                        / "sub-{subject}/anat/sub-{subject}_acq-mprage_T1w.nii.gz"
                    ),
                }
            },
            participant_filter="group == 'control'",
        )
        assert dataset.subjects == ["001"]
