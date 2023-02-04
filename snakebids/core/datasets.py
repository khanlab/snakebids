from __future__ import annotations

import itertools as it
import operator as op
from collections import UserDict
from string import Formatter
from typing import TYPE_CHECKING, Any, Iterable, Optional, Union, cast

import attr
import more_itertools as itx
from cached_property import cached_property
from typing_extensions import TypedDict

import snakebids.utils.sb_itertools as sb_it
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

    Properties
    ----------
    name
        Name of the component
    path
        Wildcard-filled path that matches the files for this component.
    zip_lists
        Dictionary where each key is a wildcard entity and each value is a list of the
        values found for that entity. Each of these lists has length equal to the number
        of images matched for this modality, so they can be zipped together to get a
        list of the wildcard values for each file.
    """

    name: str = attr.field(on_setattr=attr.setters.frozen)
    path: str = attr.field(on_setattr=attr.setters.frozen)
    zip_lists: dict[str, list[str]] = attr.field(
        on_setattr=attr.setters.frozen, converter=dict
    )

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

    @property
    def input_name(self):
        """Alias of :attr:`name <snakebids.BidsComponent.name>`

        Name of the component
        """
        return self.name

    @property
    def input_path(self):
        """Alias of :attr:`path <snakebids.BidsComponent.path>`

        Wildcard-filled path that matches the files for this component.
        """
        return self.path

    @property
    def input_zip_lists(self):
        """Alias of :attr:`zip_lists <snakebids.BidsComponent.zip_lists>`

        Dictionary where each key is a wildcard entity and each value is a list of the
        values found for that entity. Each of these lists has length equal to the number
        of images matched for this modality, so they can be zipped together to get a
        list of the wildcard values for each file.
        """
        return self.zip_lists

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
    def entities(self):
        """Component entities and their associated values

        Dictionary where each key is an entity and each value is a list of the
        unique values found for that entity. These lists might not be the same length.
        """
        if self._input_lists is None:
            self._input_lists = {
                entity: list(set(values)) for entity, values in self.zip_lists.items()
            }
        return self._input_lists

    @property_alias(entities, "entities", "snakebids.BidsComponent.entities")
    def input_lists(self):
        return self.entities

    @property
    def wildcards(self):
        """Wildcards in brace-wrapped syntax

        Dictionary where each key is the name of a wildcard entity, and each value is
        the Snakemake wildcard used for that entity.
        """
        if self._input_wildcards is None:
            self._input_wildcards = {
                entity: f"{{{entity}}}" for entity in self.zip_lists
            }
        return self._input_wildcards

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


if TYPE_CHECKING:
    _BidsComponentsType = UserDict[str, BidsComponent]
else:
    # UserDict is not subscriptable in py37
    _BidsComponentsType = UserDict


class BidsDataset(_BidsComponentsType):
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
    sessions found in the dataset
    """

    # pylint: disable=super-init-not-called
    def __init__(self, data: Any):
        self.data = dict(data)  # type: ignore

    def __setitem__(self, _: Any, __: Any):
        raise NotImplementedError(
            f"Modification of {self.__class__.__name__} is not yet supported"
        )

    @cached_property
    def path(self):
        """Dict mapping :class:`BidsComponents <snakebids.BidsComponent>` names to \
        their ``paths``.
        """
        return {key: value.path for key, value in self.data.items()}

    @property_alias(path, "path", "snakebids.BidsDataset.path")
    def input_path(self):
        return self.path

    @cached_property
    def zip_lists(self):
        """Dict mapping :class:`BidsComponents <snakebids.BidsComponent>` names to \
        their ``zip_lists``
        """
        return {key: value.zip_lists for key, value in self.data.items()}

    @property_alias(zip_lists, "zip_lists", "snakebids.BidsDataset.zip_lists")
    def input_zip_lists(self):
        return self.zip_lists

    @cached_property
    def entities(self):
        """Dict mapping :class:`BidsComponents <snakebids.BidsComponent>` names to \
        to their :attr:`entities <snakebids.BidsComponent.entities>`
        """
        return {key: value.entities for key, value in self.data.items()}

    @property_alias(entities, "entities", "snakebids.BidsDataset.entities")
    def input_lists(self):
        return self.entities

    @cached_property
    def wildcards(self):
        """Dict mapping :class:`BidsComponents <snakebids.BidsComponent>` names to \
        their :attr:`wildcards <snakebids.BidsComponent.wildcards>`
        """
        return {key: value.input_wildcards for key, value in self.data.items()}

    @property_alias(wildcards, "wildcards", "snakebids.BidsDataset.wildcards")
    def input_wildcards(self):
        return self.wildcards

    @cached_property
    def subjects(self):
        """A list of the subjects in the dataset."""
        return [
            *{
                *it.chain.from_iterable(
                    input_list["subject"]
                    for input_list in self.entities.values()
                    if "subject" in input_list
                )
            }
        ]

    @cached_property
    def sessions(self):
        """A list of the sessions in the dataset."""
        return [
            *{
                *it.chain.from_iterable(
                    input_list["session"]
                    for input_list in self.entities.values()
                    if "session" in input_list
                )
            }
        ]

    @cached_property
    def subj_wildcards(self):
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

    @property
    def as_dict(self):
        """Get the layout as a legacy dict

        Included primarily for backward compatability with older versions of snakebids,
        where generate_inputs() returned a dict rather than the `BidsDataset` class

        Returns
        -------
        BidsDatasetDict
        """
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
    def from_iterable(cls, iterable: Iterable[BidsComponent]):
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
        return cls(indexed)
