from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping


STANDARD_AA_CODES = frozenset("ACDEFGHIKLMNPQRSTVWY")
STANDARD_AA_CODES_WITH_D = STANDARD_AA_CODES | {aa.lower() for aa in STANDARD_AA_CODES}


@dataclass(frozen=True)
class UnusualResidue:
    position: int
    name: str
    before_modification: str | None = None
    note: str | None = None

    @classmethod
    def from_value(cls, value: "UnusualResidue | Mapping[str, Any]") -> "UnusualResidue":
        if isinstance(value, cls):
            return value
        modification = value.get("modificationType") or {}
        return cls(
            position=int(value["position"]),
            name=value.get("name") or modification.get("name"),
            before_modification=value.get("beforeModification"),
            note=value.get("note"),
        )

    def to_core_payload(self) -> dict[str, Any]:
        return {
            "position": self.position,
            "beforeModification": self.before_modification,
            "note": self.note,
            "modificationType": {"name": self.name},
        }


@dataclass(frozen=True)
class IntrachainBond:
    position1: int
    position2: int
    bond_type: str
    chain_participating: str = "SSB"
    cycle_type: str | None = None

    @classmethod
    def from_value(cls, value: "IntrachainBond | Mapping[str, Any]") -> "IntrachainBond":
        if isinstance(value, cls):
            return value
        bond_type = value.get("bond_type") or value.get("type")
        if isinstance(bond_type, Mapping):
            bond_type = bond_type.get("name")
        chain_participating = value.get("chain_participating") or value.get("chainParticipating") or "SSB"
        if isinstance(chain_participating, Mapping):
            chain_participating = chain_participating.get("name")
        cycle_type = value.get("cycle_type") or value.get("cycleType")
        if isinstance(cycle_type, Mapping):
            cycle_type = cycle_type.get("name")
        return cls(
            position1=int(value["position1"]),
            position2=int(value["position2"]),
            bond_type=str(bond_type),
            chain_participating=str(chain_participating),
            cycle_type=cycle_type,
        )

    def to_core_payload(self) -> dict[str, Any]:
        return {
            "position1": self.position1,
            "position2": self.position2,
            "type": {"name": self.bond_type},
            "chainParticipating": {"name": self.chain_participating},
            "cycleType": {"name": self.cycle_type} if self.cycle_type else None,
        }


@dataclass(frozen=True)
class PeptideInput:
    sequence: str
    unusual_amino_acids: tuple[UnusualResidue, ...] = field(default_factory=tuple)
    intrachain_bonds: tuple[IntrachainBond, ...] = field(default_factory=tuple)
    n_terminal: str | None = None
    c_terminal: str | None = None

    def to_api_kwargs(self) -> dict[str, Any]:
        return {
            "sequence": self.sequence,
            "unusual_amino_acids": list(self.unusual_amino_acids),
            "intrachain_bonds": list(self.intrachain_bonds),
            "n_terminal": self.n_terminal,
            "c_terminal": self.c_terminal,
        }


@dataclass(frozen=True)
class PeptideParseResult:
    sequence: str | None
    is_cyclic: bool
    cyclization: str | None
    normalized_smiles: str | None
    input_format: str
    unsupported_reason: str | None = None
