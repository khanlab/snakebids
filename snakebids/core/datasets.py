from __future__ import annotations

import itertools as it
import operator as op
import textwrap
import warnings
from math import inf
from string import Formatter
from typing import Any, Iterable, Optional, Union, cast

import attr
import more_itertools as itx
from bids import BIDSLayout
from cached_property import cached_property
from pvandyken.deprecated import deprecated
from typing_extensions import TypedDict

import snakebids.utils.sb_itertools as sb_it
from snakebids.io.console import get_console_size
from snakebids.io.printing import format_zip_lists, quote_wrap
from snakebids.types import UserDictPy37, ZipLists
from snakebids.utils.utils import property_alias


class BidsDatasetDict(TypedDict):
    """Dict equivalent of BidsInputs, for backwards-compatibility"""

    input_path: dict[str, str]
    input_zip_lists: dict[str, dict[str, list[str]]]
    input_lists: dict[str, dict[str, list[str]]]
    input_wildcards: dict[str, dict[str, str]]
    subjects: list[str]
    sessions: list[str]
    subj_wildcards: dict[str, str]


@attr.define
class BidsComponent:
    """Component of a BidsDataset mapping entities to their resolved values

    BidsComponents are immutable: their values cannot be altered.
    """

    name: str = attr.field(on_setattr=attr.setters.frozen)
    """Name of the component"""

    path: str = attr.field(on_setattr=attr.setters.frozen)
    """Wildcard-filled path that matches the files for this component."""
    zip_lists: ZipLists = attr.field(on_setattr=attr.setters.frozen, converter=dict)
    """Table of unique wildcard groupings for each member in the component.

    Dictionary where each key is a wildcard entity and each value is a list of the
    values found for that entity. Each of these lists has length equal to the number
    of images matched for this modality, so they can be zipped together to get a
    list of the wildcard values for each file.
    """

    def __repr__(self) -> str:
        return self.pformat()

    def pformat(self, max_width: int | float | None = None, tabstop: int = 4):
        width = max_width or get_console_size()[0] or inf
        body = [
            f"name={quote_wrap(self.name)},",
            f"path={quote_wrap(self.path)},",
            f"zip_lists={format_zip_lists(self.zip_lists, width - tabstop, tabstop)},",
        ]
        output = [
            "BidsComponent(",
            textwrap.indent("\n".join(body), " " * tabstop),
            ")",
        ]
        return "\n".join(output)

    @zip_lists.validator  # type: ignore
    def _validate_zip_lists(self, _, value: dict[str, list[str]]):
        lengths = {len(val) for val in value.values()}
        if len(lengths) > 1:
            raise ValueError("zip_lists must all be of equal length")
        _, raw_fields, *_ = sb_it.unpack(
            zip(*Formatter().parse(self.path)), [[], [], []]
        )
        fields = set(filter(None, raw_fields))
        if fields != set(value):
            raise ValueError(
                "zip_lists entries must match the wildcards in input_path: "
                f"{self.path}: {fields} != zip_lists: {set(value)}"
            )

    # Note: we can't use cached property here because it's incompatible with slots.
    _input_lists: Optional[dict[str, list[str]]] = attr.field(
        default=None, init=False, eq=False, repr=False
    )
    _input_wildcards: Optional[dict[str, str]] = attr.field(
        default=None, init=False, eq=False, repr=False
    )
    _entities: Optional[list[str]] = attr.field(
        default=None, init=False, eq=False, repr=False
    )

    @property
    def entities(self) -> dict[str, list[str]]:
        """Component entities and their associated values

        Dictionary where each key is an entity and each value is a list of the
        unique values found for that entity. These lists might not be the same length.
        """
        if self._input_lists is None:
            self._input_lists = {
                entity: list(set(values)) for entity, values in self.zip_lists.items()
            }
        return self._input_lists

    @property
    def wildcards(self) -> dict[str, str]:
        """Wildcards in brace-wrapped syntax

        Dictionary where each key is the name of a wildcard entity, and each value is
        the Snakemake wildcard used for that entity.
        """
        if self._input_wildcards is None:
            self._input_wildcards = {
                entity: f"{{{entity}}}" for entity in self.zip_lists
            }
        return self._input_wildcards

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

    @property
    def input_zip_lists(self) -> ZipLists:
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

    def __eq__(self, other: Union["BidsComponent", object]):
        if not isinstance(other, BidsComponent):
            return False

        if self.name != other.name:
            return False

        if self.path != other.path:
            return False

        def sorted_items(dictionary: dict[str, list[str]]):
            return sorted(dictionary.items(), key=op.itemgetter(0))

        if set(self.zip_lists) != set(other.zip_lists):
            return False

        if not other.zip_lists and not self.zip_lists:
            return True

        other_items = cast(
            "list[list[str]]", list(zip(*sorted_items(other.zip_lists)))[1]
        )
        our_items = cast("list[list[str]]", list(zip(*sorted_items(self.zip_lists)))[1])

        return set(zip(*our_items)) == set(zip(*other_items))


class BidsDataset(UserDictPy37[str, BidsComponent]):
    """A bids dataset parsed by pybids, organized into BidsComponents.

    BidsDatasets are typically generated using :func:`generate_inputs()
    <snakebids.generate_inputs>`, which reads the ``pybids_inputs`` field in your
    snakemake config file and, for each entry, creates a BidsComponent using the
    provided name, wildcards, and filters.

    Individual components can be accessed using bracket-syntax: (e.g.
    ``inputs["t1w"]``). Component access attributes along with the component name in
    brackets can also be used. For example, ``BidsDataset.entities["t1w"]`` and
    ``BidsDataset["t1w"].entities`` return the same thing.

    Provides access to summarizing information, for instance, the set of all subjects or
    sessions found in the dataset.
    """

    layout: Optional[BIDSLayout]
    """
    Underlying layout generated from pybids. Note that this will be set to None if
    custom paths are used to generate every :class:`component <BidsComponent>`
    """

    # pylint: disable=super-init-not-called
    def __init__(self, data: Any, layout: BIDSLayout | None = None):
        super().__init__(data)
        self.layout = layout

    def __setitem__(self, _: Any, __: Any):
        raise NotImplementedError(
            f"Modification of {self.__class__.__name__} is not yet supported"
        )

    def __repr__(self) -> str:
        return self.pformat()

    def pformat(self, max_width: int | float | None = None, tabstop: int = 4):
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
    def zip_lists(self) -> dict[str, ZipLists]:
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
    def entities(self) -> dict[str, dict[str, list[str]]]:
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
    def wildcards(self) -> dict[str, dict[str, str]]:
        """Dict mapping :class:`BidsComponents <snakebids.BidsComponent>` names to \
        their :attr:`wildcards <snakebids.BidsComponent.wildcards>`
        """
        return {key: value.input_wildcards for key, value in self.items()}

    @cached_property
    def subjects(self):
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
    def input_lists(self) -> dict[str, dict[str, list[str]]]:
        return self.entities

    @property_alias(zip_lists, "zip_lists", "snakebids.BidsDataset.zip_lists")
    def input_zip_lists(self) -> dict[str, dict[str, list[str]]]:
        return self.zip_lists

    @property_alias(wildcards, "wildcards", "snakebids.BidsDataset.wildcards")
    def input_wildcards(self) -> dict[str, dict[str, str]]:
        return self.wildcards

    @property
    def as_dict(self):
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
                input_lists=self.entities,
                input_wildcards=self.input_wildcards,
                input_zip_lists=self.input_zip_lists,
                subjects=self.subjects,
                sessions=self.sessions,
                subj_wildcards=self.subj_wildcards,
            )

    @classmethod
    def from_iterable(
        cls, iterable: Iterable[BidsComponent], layout: BIDSLayout | None = None
    ):
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
            raise ValueError(
                list(itx.duplicates_everseen([c.name for c in components]))
            )
        return cls(indexed, layout)
