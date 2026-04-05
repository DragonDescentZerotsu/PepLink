from __future__ import annotations

import contextlib
import io
from typing import Mapping, Sequence

from rdkit import Chem

try:
    import selfies as sf
except ImportError:  # pragma: no cover - optional at runtime
    sf = None

from ._aa_seqs_to_smiles_core import Peptide as CorePeptide

from ._chem import canonical_smiles
from .data import (
    load_registered_aa_smiles,
    load_default_c_terminal_smiles,
    load_default_n_terminal_smiles,
    merged_mapping,
)
from .exceptions import MissingDependencyError, UnsupportedPeptideError, ValidationError
from .models import IntrachainBond, PeptideInput, STANDARD_AA_CODES_WITH_D, UnusualResidue


SUPPORTED_INTRACHAIN_BONDS = frozenset(
    {
        "DSB",
        "AMD",
        "TIE",
        "DCB",
        "EST",
        "AMN",
        "p-XylB",
        "TRZB",
        "(E)-but-2-enyl-B",
        "BisMeBn-B",
        "but-2-ynyl-B",
    }
)


def _normalize_input(
    sequence: str,
    unusual_amino_acids: Sequence[UnusualResidue | Mapping[str, object]] | None,
    intrachain_bonds: Sequence[IntrachainBond | Mapping[str, object]] | None,
    n_terminal: str | None,
    c_terminal: str | None,
) -> PeptideInput:
    sequence = sequence.strip()
    if not sequence:
        raise ValidationError("sequence must not be empty")
    unusual = tuple(UnusualResidue.from_value(item) for item in (unusual_amino_acids or ()))
    bonds = tuple(IntrachainBond.from_value(item) for item in (intrachain_bonds or ()))
    return PeptideInput(
        sequence=sequence,
        unusual_amino_acids=unusual,
        intrachain_bonds=bonds,
        n_terminal=n_terminal,
        c_terminal=c_terminal,
    )


def _validate_input(
    peptide_input: PeptideInput,
    aa_smiles: Mapping[str, str],
    n_terminal_smiles: Mapping[str, str],
    c_terminal_smiles: Mapping[str, str],
) -> None:
    x_positions = [index for index, aa in enumerate(peptide_input.sequence, start=1) if aa in {"X", "x"}]
    if len(x_positions) != len(peptide_input.unusual_amino_acids):
        raise ValidationError(
            f"sequence contains {len(x_positions)} X/x placeholders but {len(peptide_input.unusual_amino_acids)} unusual_amino_acids were provided."
        )
    unusual_positions = sorted(item.position for item in peptide_input.unusual_amino_acids)
    if unusual_positions != x_positions:
        raise ValidationError(
            f"unusual_amino_acids positions {unusual_positions} do not match X/x placeholder positions {x_positions}."
        )
    unknown_residues = sorted({aa for aa in peptide_input.sequence if aa not in {"X", "x"} and aa not in aa_smiles})
    if unknown_residues:
        raise UnsupportedPeptideError(f"Unsupported residue codes in sequence: {', '.join(unknown_residues)}")
    missing_unusual = sorted({item.name for item in peptide_input.unusual_amino_acids if item.name not in aa_smiles})
    if missing_unusual:
        raise UnsupportedPeptideError(
            "Unsupported non-canonical residues: "
            + ", ".join(missing_unusual)
            + ". Add them through aa_overrides or extend the bundled mapping."
        )
    unsupported_bonds = sorted({bond.bond_type for bond in peptide_input.intrachain_bonds if bond.bond_type not in SUPPORTED_INTRACHAIN_BONDS})
    if unsupported_bonds:
        raise UnsupportedPeptideError(
            "Unsupported intrachain bond types for PepLink v1: " + ", ".join(unsupported_bonds)
        )
    if peptide_input.n_terminal and peptide_input.n_terminal not in n_terminal_smiles:
        raise UnsupportedPeptideError(
            f"Unsupported N-terminal modification {peptide_input.n_terminal!r}. Add it through n_terminal_overrides."
        )
    if peptide_input.c_terminal and peptide_input.c_terminal not in c_terminal_smiles:
        raise UnsupportedPeptideError(
            f"Unsupported C-terminal modification {peptide_input.c_terminal!r}. Add it through c_terminal_overrides."
        )


def _build_core_peptide(
    peptide_input: PeptideInput,
    aa_smiles: Mapping[str, str],
    n_terminal_smiles: Mapping[str, str],
    c_terminal_smiles: Mapping[str, str],
) -> CorePeptide:
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        peptide = CorePeptide(
            peptide_input.sequence,
            aa_smiles_dict=dict(aa_smiles),
            intrachain_bonds=[bond.to_core_payload() for bond in peptide_input.intrachain_bonds],
            interchain_bonds=[],
            unusual_aas=[item.to_core_payload() for item in peptide_input.unusual_amino_acids],
            cTerminus=peptide_input.c_terminal,
            nTerminus=peptide_input.n_terminal,
            c_terminus_modify_name_smiles=dict(c_terminal_smiles),
            n_terminus_modify_name_smiles=dict(n_terminal_smiles),
        )
    if peptide.noise_data_flag:
        message = "; ".join(line.strip() for line in buffer.getvalue().splitlines() if line.strip())
        raise UnsupportedPeptideError(message or "The peptide builder reported an unsupported peptide definition.")
    if not peptide.ncTerminus_modified_mols:
        raise ValidationError("The peptide builder did not produce an output molecule.")
    return peptide


def peptide_to_mol(
    sequence: str,
    *,
    unusual_amino_acids: Sequence[UnusualResidue | Mapping[str, object]] | None = None,
    intrachain_bonds: Sequence[IntrachainBond | Mapping[str, object]] | None = None,
    n_terminal: str | None = None,
    c_terminal: str | None = None,
    aa_overrides: Mapping[str, str] | None = None,
    n_terminal_overrides: Mapping[str, str] | None = None,
    c_terminal_overrides: Mapping[str, str] | None = None,
) -> Chem.Mol:
    peptide_input = _normalize_input(sequence, unusual_amino_acids, intrachain_bonds, n_terminal, c_terminal)
    aa_smiles = merged_mapping(load_registered_aa_smiles(), aa_overrides)
    n_terminal_smiles = merged_mapping(load_default_n_terminal_smiles(), n_terminal_overrides)
    c_terminal_smiles = merged_mapping(load_default_c_terminal_smiles(), c_terminal_overrides)
    _validate_input(peptide_input, aa_smiles, n_terminal_smiles, c_terminal_smiles)
    peptide = _build_core_peptide(peptide_input, aa_smiles, n_terminal_smiles, c_terminal_smiles)
    return Chem.Mol(peptide.ncTerminus_modified_mols[0])


def aa_seqs_to_smiles(
    sequence: str,
    *,
    unusual_amino_acids: Sequence[UnusualResidue | Mapping[str, object]] | None = None,
    intrachain_bonds: Sequence[IntrachainBond | Mapping[str, object]] | None = None,
    n_terminal: str | None = None,
    c_terminal: str | None = None,
    output_format: str = "smiles",
    aa_overrides: Mapping[str, str] | None = None,
    n_terminal_overrides: Mapping[str, str] | None = None,
    c_terminal_overrides: Mapping[str, str] | None = None,
) -> str:
    mol = peptide_to_mol(
        sequence,
        unusual_amino_acids=unusual_amino_acids,
        intrachain_bonds=intrachain_bonds,
        n_terminal=n_terminal,
        c_terminal=c_terminal,
        aa_overrides=aa_overrides,
        n_terminal_overrides=n_terminal_overrides,
        c_terminal_overrides=c_terminal_overrides,
    )
    smiles = canonical_smiles(mol)
    if output_format == "smiles":
        return smiles
    if output_format == "selfies":
        if sf is None:
            raise MissingDependencyError("SELFIES output requires the optional 'selfies' dependency.")
        return sf.encoder(smiles)
    raise ValidationError("output_format must be either 'smiles' or 'selfies'")


def head_to_tail_bond(sequence: str) -> list[IntrachainBond]:
    if not sequence:
        raise ValidationError("sequence must not be empty")
    return [IntrachainBond(position1=1, position2=len(sequence), bond_type="AMD", chain_participating="MMB")]
