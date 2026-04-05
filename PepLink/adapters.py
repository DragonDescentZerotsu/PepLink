from __future__ import annotations

from typing import Any, Mapping

from .exceptions import UnsupportedPeptideError, ValidationError
from .models import IntrachainBond, PeptideInput, UnusualResidue


def from_dbaasp_record(record: Mapping[str, Any]) -> PeptideInput:
    record_id = record.get("id", "<unknown>")
    complexity = (record.get("complexity") or {}).get("name")
    if complexity != "Monomer":
        raise UnsupportedPeptideError(
            f"DBAASP record {record_id} has complexity {complexity!r}; PepLink v1 only supports monomer peptides."
        )
    if record.get("interchainBonds"):
        raise UnsupportedPeptideError(f"DBAASP record {record_id} contains interchain bonds, which are not supported in PepLink v1.")
    if record.get("coordinationBonds"):
        raise UnsupportedPeptideError(
            f"DBAASP record {record_id} contains coordination bonds, which are not supported in PepLink v1."
        )
    sequence = (record.get("sequence") or "").strip()
    if not sequence:
        raise ValidationError(f"DBAASP record {record_id} does not contain a monomer sequence.")
    unusual = tuple(UnusualResidue.from_value(item) for item in (record.get("unusualAminoAcids") or []))
    bonds = tuple(IntrachainBond.from_value(item) for item in (record.get("intrachainBonds") or []))
    return PeptideInput(
        sequence=sequence,
        unusual_amino_acids=unusual,
        intrachain_bonds=bonds,
        n_terminal=(record.get("nTerminus") or {}).get("name"),
        c_terminal=(record.get("cTerminus") or {}).get("name"),
    )
