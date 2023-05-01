from __future__ import annotations

import itertools as it
import textwrap
import warnings
from math import inf
from pathlib import Path
from string import Formatter
from typing import Any, Iterable, NoReturn, Optional, Sequence, overload

import attr
import more_itertools as itx
from bids import BIDSLayout
from cached_property import cached_property
from pvandyken.deprecated import deprecated
from snakemake.io import expand as sn_expand
from typing_extensions import Self, TypedDict

import snakebids.utils.sb_itertools as sb_it
from snakebids.core.filtering import filter_list
from snakebids.exceptions import DuplicateComponentError
from snakebids.io.console import get_console_size
from snakebids.io.printing import format_zip_lists, quote_wrap
from snakebids.types import UserDictPy37, ZipList
from snakebids.utils.utils import (
    ImmutableList,
    MultiSelectDict,
    get_wildcard_dict,
    property_alias,
    zip_list_eq,
)


class BidsDatasetDict(TypedDict):
    """Dict equivalent of BidsInputs, for backwards-compatibility"""

    input_path: dict[str, str]
    input_zip_lists: dict[str, dict[str, list[str]]]
    input_lists: dict[str, dict[str, list[str]]]
    input_wildcards: dict[str, dict[str, str]]
    subjects: list[str]
    sessions: list[str]
    subj_wildcards: dict[str, str]


class BidsComponentRow(ImmutableList[str]):
    def __init__(self, __iterable: Iterable[str], entity: str):
        super().__init__(__iterable)
        self.entity = entity

    @property
    def entities(self) -> tuple[str]:
        """The unique values associated with the component"""
        return tuple(set(self._data))

    @property
    def wildcards(self) -> str:
        """The entity name wrapped in wildcard braces"""
        return f"{{{self.entity}}}"

    @property
    def zip_lists(self) -> ZipList:
        """
        Dictionary where each key is a wildcard entity and each value is a list of the
        values found for that entity. Each of these lists has length equal to the number
        of images matched for this modality, so they can be zipped together to get a
        list of the wildcard values for each file.
        """
        return MultiSelectDict({self.entity: list(self._data)})

    def __eq__(self, other: BidsComponentRow | object) -> bool:
        if not isinstance(other, self.__class__):
            return False

        return super().__eq__(other)

    def expand(
        self,
        paths: Iterable[Path | str] | Path | str,
        allow_missing: bool = False,
        **wildcards: str | Iterable[str],
    ) -> list[str]:
        """Safely expand over given paths with component wildcards

        Uses the entity-value combinations found in the dataset to expand over the given
        paths. If no path is provided, expands over the component
        :attr:`~snakebids.BidsComponent.path` (thus returning the original files used to
        create the component). Extra wildcards can be specifed as keyword arguments.

        By default, expansion over paths with extra wildcards not accounted for by the
        component causes an error. This prevents accidental partial expansion. To allow
        the passage of extra wildcards without expansion,set ``allow_missing`` to
        ``True``.

        Uses the snakemake :ref:`expand <snakemake:snakefiles_expand>` under the hood.
        """
        return sn_expand(
            list(itx.always_iterable(paths)),
            allow_missing=allow_missing,
            **{self.entity: self._data},
            **{
                wildcard: list(itx.always_iterable(v))
                for wildcard, v in wildcards.items()
            },
        )

    def filter(
        self, *, regex_search: bool = False, **filters: str | Sequence[str]
    ) -> Self:
        """Filter component based on provided entity filters

        This method allows you to expand over a subset of your wildcards. This could be
        useful for extracting subjects from a specific patient group, running different
        rules on different aquisitions, and any other reason you may need to filter your
        data after the workflow has already started.

        Takes entities as keyword arguments assigned to values or list of values to
        select from the component. Only columns containing the provided entity-values
        are kept. If no matches are found, a component with the all the original
        entities but with no values will be returned.

        Returns a brand new :class:`~snakebids.BidsComponent`. The original component is
        not modified.

        Parameters
        ----------
        regex_search
            Treat filters as regex patterns when matching with entity-values.
        filters
            Each keyword should be the name of an entity in the component. Entities not
            found in the component will be ignored. Keywords take values or a list of
            values to be matched with the component
            :attr:`~snakebids.BidsComponent.zip_lists`
        """
        entity, data = itx.first(
            filter_list(
                {self.entity: self._data}, filters, regex_search=regex_search
            ).items()
        )
        return self.__class__(data, entity=entity)


@attr.define(kw_only=True)
class BidsPartialComponent:
    """Primitive representation of a bids data component

    See :class:`BidsComponent` for an extended definition of a data component.

    ``BidsPartialComponents`` are typically derived from :class:`BidsComponent`. They do
    not store path information, and do not represent *real* data, only a table of
    entity-values.

    Despite this, ``BidsPartialComponents`` still allow you to expand the data table
    over new paths, allowing you to derive paths from your source dataset.

    ``BidsPartialComponents`` are immutable: their values cannot be altered.
    """

    zip_lists: ZipList = attr.field(
        on_setattr=attr.setters.frozen, converter=MultiSelectDict
    )
    """Table of unique wildcard groupings for each member in the component.

    Dictionary where each key is a wildcard entity and each value is a list of the
    values found for that entity. Each of these lists has length equal to the number
    of images matched for this modality, so they can be zipped together to get a
    list of the wildcard values for each file.
    """

    def __repr__(self) -> str:
        return self.pformat()

    @overload
    def __getitem__(self, __key: str) -> BidsComponentRow:
        ...

    @overload
    def __getitem__(self, __key: tuple[str, ...]) -> BidsPartialComponent:
        ...

    def __getitem__(
        self, __key: str | tuple[str, ...]
    ) -> BidsComponentRow | BidsPartialComponent:
        if isinstance(__key, tuple):
            # Use dict.fromkeys for de-duplication
            return BidsPartialComponent(
                zip_lists={key: self.zip_lists[key] for key in dict.fromkeys(__key)}
            )
        return BidsComponentRow(self.zip_lists[__key], entity=__key)

    def __bool__(self) -> bool:
        """Truth of a BidsComponent is based on whether it has values

        It is not based on whether it has any entities. This is because
        :meth:`~BidsPartialComponent.filter` returns a component retaining all entities
        but with no values, making that the standard for emptiness. It also makes it
        consistent with :class:`BidsComponentRow`, which always has an entity name
        stored, but may or may not have values.
        """
        return bool(next(iter(self.zip_lists)))

    def _pformat_body(self) -> None | str | list[str]:
        return None

    def pformat(self, max_width: int | float | None = None, tabstop: int = 4) -> str:
        width = max_width or get_console_size()[0] or inf
        body = it.chain(
            itx.always_iterable(self._pformat_body() or []),
            [
                "zip_lists="
                f"{format_zip_lists(self.zip_lists, width - tabstop, tabstop)},",
            ],
        )
        output = [
            f"{self.__class__.__name__}(",
            textwrap.indent("\n".join(body), " " * tabstop),
            ")",
        ]
        return "\n".join(output)

    @zip_lists.validator  # type: ignore
    def _validate_zip_lists(self, __attr: str, value: dict[str, list[str]]) -> None:
        lengths = {len(val) for val in value.values()}
        if len(lengths) > 1:
            raise ValueError("zip_lists must all be of equal length")

    # Note: we can't use cached property here because it's incompatible with slots.
    _input_lists: Optional[MultiSelectDict[str, list[str]]] = attr.field(
        default=None, init=False, eq=False, repr=False
    )
    _input_wildcards: Optional[MultiSelectDict[str, str]] = attr.field(
        default=None, init=False, eq=False, repr=False
    )
    _entities: Optional[list[str]] = attr.field(
        default=None, init=False, eq=False, repr=False
    )

    @property
    def entities(self) -> MultiSelectDict[str, list[str]]:
        """Component entities and their associated values

        Dictionary where each key is an entity and each value is a list of the
        unique values found for that entity. These lists might not be the same length.
        """
        if self._input_lists is None:
            self._input_lists = MultiSelectDict(
                {entity: list(set(values)) for entity, values in self.zip_lists.items()}
            )
        return self._input_lists

    @property
    def wildcards(self) -> MultiSelectDict[str, str]:
        """Wildcards in brace-wrapped syntax

        Dictionary where each key is the name of a wildcard entity, and each value is
        the Snakemake wildcard used for that entity.
        """
        if self._input_wildcards is None:
            self._input_wildcards = MultiSelectDict(get_wildcard_dict(self.zip_lists))
        return self._input_wildcards

    @property
    def input_zip_lists(self) -> ZipList:
        """Alias of :attr:`zip_lists <snakebids.BidsComponent.zip_lists>`

        Dictionary where each key is a wildcard entity and each value is a list of the
        values found for that entity. Each of these lists has length equal to the number
        of images matched for this modality, so they can be zipped together to get a
        list of the wildcard values for each file.
        """
        return self.zip_lists

    @property_alias(entities, "entities", "snakebids.BidsComponent.entities")
    def input_lists(self):
        return self.entities

    @property_alias(wildcards, "wildcards", "snakebids.BidsComponent.wildcards")
    def input_wildcards(self):
        return self.wildcards

    def __eq__(self, other: BidsComponent | object) -> bool:
        if not isinstance(other, self.__class__):
            return False

        return zip_list_eq(self.zip_lists, other.zip_lists)

    def expand(
        self,
        paths: Iterable[Path | str] | Path | str,
        allow_missing: bool = False,
        **wildcards: str | Iterable[str],
    ) -> list[str]:
        """Safely expand over given paths with component wildcards

        Uses the entity-value combinations found in the dataset to expand over the given
        paths. If no path is provided, expands over the component
        :attr:`~snakebids.BidsComponent.path` (thus returning the original files used to
        create the component). Extra wildcards can be specifed as keyword arguments.

        By default, expansion over paths with extra wildcards not accounted for by the
        component causes an error. This prevents accidental partial expansion. To allow
        the passage of extra wildcards without expansion,set ``allow_missing`` to
        ``True``.

        Uses the snakemake :ref:`expand <snakemake:snakefiles_expand>` under the hood.
        """
        inner_expand = sn_expand(
            list(itx.always_iterable(paths)),
            zip,
            allow_missing=True if wildcards else allow_missing,
            **self.zip_lists,
        )
        if not wildcards:
            return inner_expand

        return sn_expand(
            inner_expand,
            allow_missing=allow_missing,
            # Turn all the wildcard items into lists because Snakemake doesn't handle
            # iterables very well
            **{
                wildcard: list(itx.always_iterable(v))
                for wildcard, v in wildcards.items()
            },
        )

    def filter(
        self, *, regex_search: bool = False, **filters: str | Sequence[str]
    ) -> Self:
        """Filter component based on provided entity filters

        This method allows you to expand over a subset of your wildcards. This could be
        useful for extracting subjects from a specific patient group, running different
        rules on different aquisitions, and any other reason you may need to filter your
        data after the workflow has already started.

        Takes entities as keyword arguments assigned to values or list of values to
        select from the component. Only columns containing the provided entity-values
        are kept. If no matches are found, a component with the all the original
        entities but with no values will be returned.

        Returns a brand new :class:`~snakebids.BidsComponent`. The original component is
        not modified.

        Parameters
        ----------
        regex_search
            Treat filters as regex patterns when matching with entity-values.
        filters
            Each keyword should be the name of an entity in the component. Entities not
            found in the component will be ignored. Keywords take values or a list of
            values to be matched with the component
            :attr:`~snakebids.BidsComponent.zip_lists`
        """
        return attr.evolve(
            self,
            zip_lists=filter_list(self.zip_lists, filters, regex_search=regex_search),
        )


@attr.define(kw_only=True)
class BidsComponent(BidsPartialComponent):
    """Representation of a bids data component

    A component is a set of data entries all corresponding to the same type of object.
    Entries vary over a set of entities. For example, a component may represent all the
    unprocessed, T1-weighted anatomical images aqcuired from a group of 100 subjects,
    across 2 sessions, with three runs per session. Here, the subject, session, and run
    are the entities over which the component varies. Each entry in the component has
    a single value assigned for each of the three entities (e.g subject 002, session
    01, run 1).

    Each entry can be defined soley by its wildcard values. The complete collection of
    entries can thus be stored as a table, where each row represents an entity and each
    column represents an entry.

    ``BidsComponent`` stores and indexes this table. It uses 'row-first' indexing,
    meaning first an entity is selected, then an entry. It also has a number of
    properties and methods making it easier to incorporate the data in a snakemake
    workflow.

    In addition, ``BidsComponent`` stores a template :attr:`~BidsComponent.path <path>`
    derived from the source dataset. This path is used by the
    :meth:`~BidsComponent.expand` method to recreate the original filesystem paths.

    The real power of the ``BidsComponent``, however, is in creating derived paths based
    on the original dataset. Using the :meth`~BidsComponent.expand` method, you can pass
    new paths with ``{wildcard}`` placeholders wrapped in braces and named according to
    the entities in the component. These placeholders will be substituted with the
    entity values saved in the table, giving you a list of paths the same length as the
    number of entries in the component.

    BidsComponents are immutable: their values cannot be altered.
    """

    name: str = attr.field(on_setattr=attr.setters.frozen)
    """Name of the component"""

    path: str = attr.field(on_setattr=attr.setters.frozen)
    """Wildcard-filled path that matches the files for this component."""

    zip_lists: ZipList = attr.field(
        on_setattr=attr.setters.frozen, converter=MultiSelectDict
    )
    """Table of unique wildcard groupings for each member in the component.

    Dictionary where each key is a wildcard entity and each value is a list of the
    values found for that entity. Each of these lists has length equal to the number
    of images matched for this modality, so they can be zipped together to get a
    list of the wildcard values for each file.
    """

    def __repr__(self) -> str:
        return self.pformat()

    def _pformat_body(self):
        return [
            f"name={quote_wrap(self.name)},",
            f"path={quote_wrap(self.path)},",
        ]

    @zip_lists.validator  # type: ignore
    def _validate_zip_lists(self, __attr: str, value: dict[str, list[str]]) -> None:
        super()._validate_zip_lists(__attr, value)
        _, raw_fields, *_ = sb_it.unpack(
            zip(*Formatter().parse(self.path)), [[], [], []]
        )
        fields = set(filter(None, raw_fields))
        if fields != set(value):
            raise ValueError(
                "zip_lists entries must match the wildcards in input_path: "
                f"{self.path}: {fields} != zip_lists: {set(value)}"
            )

    @property
    def input_name(self) -> str:
        """Alias of :attr:`name <snakebids.BidsComponent.name>`

        Name of the component
        """
        return self.name

    @property
    def input_path(self) -> str:
        """Alias of :attr:`path <snakebids.BidsComponent.path>`

        Wildcard-filled path that matches the files for this component.
        """
        return self.path

    def __eq__(self, other: BidsComponent | object) -> bool:
        if not isinstance(other, self.__class__):
            return False

        if self.name != other.name:
            return False

        if self.path != other.path:
            return False

        return super().__eq__(other)

    def expand(
        self,
        paths: Iterable[Path | str] | Path | str | None = None,
        allow_missing: bool = False,
        **wildcards: str | Iterable[str],
    ) -> list[str]:
        """Safely expand over given paths with component wildcards

        Uses the entity-value combinations found in the dataset to expand over the given
        paths. If no path is provided, expands over the component
        :attr:`~snakebids.BidsComponent.path` (thus returning the original files used to
        create the component). Extra wildcards can be specifed as keyword arguments.

        By default, expansion over paths with extra wildcards not accounted for by the
        component causes an error. This prevents accidental partial expansion. To allow
        the passage of extra wildcards without expansion,set ``allow_missing`` to
        ``True``.

        Uses the snakemake :ref:`expand <snakemake:snakefiles_expand>` under the hood.
        """
        paths = paths or self.path
        return super().expand(paths, allow_missing, **wildcards)


class BidsDataset(UserDictPy37[str, BidsComponent]):
    """A bids dataset parsed by pybids, organized into BidsComponents.

    BidsDatasets are typically generated using :func:`generate_inputs()
    <snakebids.generate_inputs>`, which reads the ``pybids_inputs`` field in your
    snakemake config file and, for each entry, creates a BidsComponent using the
    provided name, wildcards, and filters.

    Individual components can be accessed using bracket-syntax: (e.g.
    ``inputs["t1w"]``).

    Provides access to summarizing information, for instance, the set of all subjects or
    sessions found in the dataset.
    """

    layout: Optional[BIDSLayout]
    """
    Underlying layout generated from pybids. Note that this will be set to None if
    custom paths are used to generate every :class:`component <BidsComponent>`
    """

    def __init__(self, data: Any, layout: BIDSLayout | None = None) -> None:
        super().__init__(data)
        self.layout = layout

    def __getitem__(self, key: str) -> BidsComponent:
        try:
            return super().__getitem__(key)
        except KeyError as err:
            if key in {
                "input_path",
                "input_zip_lists",
                "input_lists",
                "input_wildcards",
                "subjects",
                "sessions",
                "subj_wildcards",
            }:
                raise KeyError(
                    "As of v0.8, generate_inputs() no longer returns a dict by "
                    f"default, but an instance of BidsDataset. As such, '{key}' can no "
                    "longer be accessed via brackets '[]' as before. The original dict "
                    "can be returned by setting `use_bids_inputs` to False in the call "
                    "to generate_inputs(). However, we encourage you to transition to "
                    "the use of `BidsDataset` for long term support"
                ) from err
            raise err

    def __setitem__(self, _: Any, __: Any) -> NoReturn:
        raise NotImplementedError(
            f"Modification of {self.__class__.__name__} is not yet supported"
        )

    def __repr__(self) -> str:
        return self.pformat()

    def pformat(self, max_width: int | float | None = None, tabstop: int = 4) -> str:
        width = max_width or get_console_size()[0] or inf
        body = [
            f"{quote_wrap(name)}: {comp.pformat(width - tabstop, tabstop)},"
            for name, comp in self.items()
        ]
        output = [
            "BidsDataset({",
            textwrap.indent("\n".join(body), " " * tabstop),
            "})",
        ]
        return "\n".join(output)

    @cached_property
    @deprecated(
        details="""
        The behaviour of path will change in an upcoming release, where it will refer
        instead to the root path of the dataset. Please access component paths using
        :attr:`Dataset[\\<component_name\\>].path <BidsComponent.path>`
        """,
        deprecated_in="0.8.0",
        admonition="warning",
    )
    def path(self) -> dict[str, str]:
        """Dict mapping :class:`BidsComponents <snakebids.BidsComponent>` names to \
        their ``paths``.
        """
        return {key: value.path for key, value in self.items()}

    @cached_property
    @deprecated(
        details="""
        The behaviour of zip_lists will change in an upcoming release, where it will
        refer instead to the consensus of entity groups across all components in the
        dataset. Please access component zip_lists using
        :attr:`Dataset[\\<component_name\\>].zip_lists <BidsComponent.zip_lists>`
        """,
        deprecated_in="0.8.0",
        admonition="warning",
    )
    def zip_lists(self) -> dict[str, ZipList]:
        """Dict mapping :class:`BidsComponents <snakebids.BidsComponent>` names to \
        their ``zip_lists``
        """
        return {key: value.zip_lists for key, value in self.items()}

    @cached_property
    @deprecated(
        details="""
        The behaviour of entities will change in the 1.0 release, where it will refer
        instead to the union of all entity-values across all components in the dataset.
        Please access component entity lists using
        :attr:`Dataset[\\<component_name\\>].entities <BidsComponent.entities>`
        """,
        deprecated_in="0.8.0",
        admonition="warning",
    )
    def entities(self) -> dict[str, MultiSelectDict[str, list[str]]]:
        """Dict mapping :class:`BidsComponents <snakebids.BidsComponent>` names to \
        to their :attr:`entities <snakebids.BidsComponent.entities>`
        """
        return {key: value.entities for key, value in self.items()}

    @cached_property
    @deprecated(
        details="""
        The behaviour of wildcards will change in an upcoming release, where it will
        refer instead to the union of all entity-wildcard mappings across all components
        in the dataset. Please access component wildcards using
        :attr:`Dataset[\\<component_name\\>].wildcards <BidsComponent.wildcards>`
        """,
        deprecated_in="0.8.0",
        admonition="warning",
    )
    def wildcards(self) -> dict[str, MultiSelectDict[str, str]]:
        """Dict mapping :class:`BidsComponents <snakebids.BidsComponent>` names to \
        their :attr:`wildcards <snakebids.BidsComponent.wildcards>`
        """
        return {key: value.input_wildcards for key, value in self.items()}

    @cached_property
    def subjects(self) -> list[str]:
        """A list of the subjects in the dataset."""
        return [
            *{
                *it.chain.from_iterable(
                    component.entities.get("subject", []) for component in self.values()
                )
            }
        ]

    @cached_property
    def sessions(self) -> list[str]:
        """A list of the sessions in the dataset."""
        return [
            *{
                *it.chain.from_iterable(
                    component.entities.get("session", []) for component in self.values()
                )
            }
        ]

    @cached_property
    def subj_wildcards(self) -> dict[str, str]:
        """The subject and session wildcards applicable to this dataset.

        ``{"subject":"{subject}"}`` if there is only one session, ``{"subject":
        "{subject}", "session": "{session}"}`` if there are multiple sessions.
        """
        if len(self.sessions) == 0:
            return {"subject": "{subject}"}
        return {
            "subject": "{subject}",
            "session": "{session}",
        }

    @property_alias(path, "path", "snakebids.BidsDataset.path")
    def input_path(self) -> dict[str, str]:
        return self.path

    @property_alias(entities, "entities", "snakebids.BidsDataset.entities")
    def input_lists(self) -> dict[str, MultiSelectDict[str, list[str]]]:
        return self.entities

    @property_alias(zip_lists, "zip_lists", "snakebids.BidsDataset.zip_lists")
    def input_zip_lists(self) -> dict[str, MultiSelectDict[str, list[str]]]:
        return self.zip_lists

    @property_alias(wildcards, "wildcards", "snakebids.BidsDataset.wildcards")
    def input_wildcards(self) -> dict[str, MultiSelectDict[str, str]]:
        return self.wildcards

    @property
    def as_dict(self) -> BidsDatasetDict:
        """Get the layout as a legacy dict

        Included primarily for backward compatability with older versions of snakebids,
        where generate_inputs() returned a dict rather than the `BidsDataset` class

        Returns
        -------
        BidsDatasetDict
        """
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            return BidsDatasetDict(
                input_path=self.input_path,
                input_lists={k: dict(val) for k, val in self.entities.items()},
                input_wildcards={
                    k: dict(val) for k, val in self.input_wildcards.items()
                },
                input_zip_lists={
                    k: dict(val) for k, val in self.input_zip_lists.items()
                },
                subjects=self.subjects,
                sessions=self.sessions,
                subj_wildcards=self.subj_wildcards,
            )

    @classmethod
    def from_iterable(
        cls, iterable: Iterable[BidsComponent], layout: BIDSLayout | None = None
    ) -> BidsDataset:
        """Construct Dataset from iterable of BidsComponents

        Parameters
        ----------
        iterable : Iterable[BidsComponent]

        Returns
        -------
        BidsDataset
        """
        components = list(iterable)
        indexed = {bidsinput.name: bidsinput for bidsinput in components}
        if not len(components) == len(indexed):
            raise DuplicateComponentError(
                list(itx.duplicates_everseen([c.name for c in components]))
            )
        return cls(indexed, layout)
