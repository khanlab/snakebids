"""Reading and querying the ``participants.tsv`` file of a BIDS dataset."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from snakebids.exceptions import ConfigError

_logger = logging.getLogger(__name__)

PARTICIPANTS_FILE = "participants.tsv"


def read_participants_tsv(bids_dir: Path | str) -> pd.DataFrame | None:
    """Read the ``participants.tsv`` file at the root of a BIDS dataset.

    The table is indexed by subject label (``participant_id`` with the ``sub-`` prefix
    removed), and the index is named ``subject``. ``n/a`` values are parsed as missing
    values, and the types of all other columns are inferred.

    Parameters
    ----------
    bids_dir
        Root of the BIDS dataset

    Returns
    -------
    pd.DataFrame | None
        The participants table, or None if the dataset has no ``participants.tsv``

    Raises
    ------
    ConfigError
        If the file cannot be parsed or has no ``participant_id`` column
    """
    path = Path(bids_dir, PARTICIPANTS_FILE)
    if not path.is_file():
        return None
    try:
        table = pd.read_csv(
            path,
            sep="\t",
            dtype={"participant_id": str},
            na_values=["n/a"],
            keep_default_na=False,
        )
    except (pd.errors.ParserError, pd.errors.EmptyDataError) as err:
        msg = f"Could not parse {path}: {err}"
        raise ConfigError(msg) from err
    if "participant_id" not in table.columns:
        msg = f"{path} does not have the required 'participant_id' column"
        raise ConfigError(msg)
    table.index = pd.Index(
        table["participant_id"].str.removeprefix("sub-"), name="subject"
    )
    if table.index.has_duplicates:
        _logger.warning(
            "%s has duplicate participant_id entries: %s",
            path,
            ", ".join(table.index[table.index.duplicated()].unique()),
        )
    return table


def query_participants(table: pd.DataFrame, query: str) -> list[str]:
    """Get the labels of the subjects matching a query on the participants table.

    Parameters
    ----------
    table
        Participants table, as returned by :func:`read_participants_tsv`
    query
        Query string, evaluated with :meth:`pandas.DataFrame.query`, e.g.
        ``"age > 30 and group == 'patient'"``. The subject label can be referenced as
        ``subject``. Column names with spaces or special characters must be
        surrounded by backticks.

    Returns
    -------
    list[str]
        Labels (without ``sub-``) of the matching subjects

    Raises
    ------
    ConfigError
        If the query cannot be evaluated
    """
    try:
        selected = table.query(query)
    except Exception as err:
        msg = (
            f"Could not evaluate participant filter {query!r}: {err}\n"
            f"Available columns: {', '.join(map(str, table.columns))}"
        )
        raise ConfigError(msg) from err
    return list(dict.fromkeys(selected.index))
