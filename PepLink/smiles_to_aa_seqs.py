from __future__ import annotations

from rdkit import Chem

try:
    import selfies as sf
except ImportError:  # pragma: no cover - optional at runtime
    sf = None

from ._smiles_to_aa_seqs_core import smiles_to_pepseq as core_smiles_to_pepseq

from ._chem import canonical_smiles
from .exceptions import MissingDependencyError, PepLinkError, ValidationError
from .aa_seqs_to_smiles import aa_seqs_to_smiles, head_to_tail_bond
from .models import PeptideParseResult, STANDARD_AA_CODES_WITH_D


ARG_TAUTOMER_TRANSLATION = str.maketrans({"B": "R", "b": "r"})


def _canonical_cyclic_sequence(sequence: str) -> str:
    rotations = [sequence[index:] + sequence[:index] for index in range(len(sequence))]
    return min(rotations)


def _decode_input(text: str, input_format: str) -> tuple[str | None, str]:
    if input_format not in {"auto", "smiles", "selfies"}:
        raise ValidationError("input_format must be 'auto', 'smiles', or 'selfies'")
    if input_format == "smiles":
        mol = Chem.MolFromSmiles(text)
        return (canonical_smiles(mol), "smiles") if mol is not None else (None, "smiles")
    if input_format == "selfies":
        if sf is None:
            raise MissingDependencyError("SELFIES parsing requires the optional 'selfies' dependency.")
        try:
            decoded = sf.decoder(text)
        except Exception:
            return None, "selfies"
        mol = Chem.MolFromSmiles(decoded)
        return (canonical_smiles(mol), "selfies") if mol is not None else (None, "selfies")

    mol = Chem.MolFromSmiles(text)
    if mol is not None:
        return canonical_smiles(mol), "smiles"
    if sf is None:
        return None, "auto"
    try:
        decoded = sf.decoder(text)
    except Exception:
        return None, "auto"
    mol = Chem.MolFromSmiles(decoded)
    return (canonical_smiles(mol), "selfies") if mol is not None else (None, "auto")


def smiles_to_aa_seqs(text: str, *, input_format: str = "auto") -> PeptideParseResult:
    normalized_smiles, resolved_format = _decode_input(text, input_format)
    if normalized_smiles is None:
        return PeptideParseResult(
            sequence=None,
            is_cyclic=False,
            cyclization=None,
            normalized_smiles=None,
            input_format=resolved_format,
            unsupported_reason="Input could not be parsed as a valid SMILES or SELFIES string.",
        )

    _, candidate = core_smiles_to_pepseq(normalized_smiles)
    if candidate is None:
        return PeptideParseResult(
            sequence=None,
            is_cyclic=False,
            cyclization=None,
            normalized_smiles=normalized_smiles,
            input_format=resolved_format,
            unsupported_reason="PepLink v1 only reverse-parses standard amino-acid peptides with linear or head-to-tail backbones.",
        )

    is_cyclic = candidate.startswith("cyclo-")
    sequence = candidate.removeprefix("cyclo-").translate(ARG_TAUTOMER_TRANSLATION)
    if is_cyclic:
        sequence = _canonical_cyclic_sequence(sequence)
    if not sequence or any(aa not in STANDARD_AA_CODES_WITH_D for aa in sequence):
        return PeptideParseResult(
            sequence=None,
            is_cyclic=False,
            cyclization=None,
            normalized_smiles=normalized_smiles,
            input_format=resolved_format,
            unsupported_reason="Reverse parsing detected non-standard or unsupported residues.",
        )

    try:
        reconstructed = aa_seqs_to_smiles(
            sequence,
            intrachain_bonds=head_to_tail_bond(sequence) if is_cyclic else None,
        )
    except PepLinkError as exc:
        return PeptideParseResult(
            sequence=None,
            is_cyclic=False,
            cyclization=None,
            normalized_smiles=normalized_smiles,
            input_format=resolved_format,
            unsupported_reason=str(exc),
        )

    if reconstructed != normalized_smiles:
        return PeptideParseResult(
            sequence=None,
            is_cyclic=False,
            cyclization=None,
            normalized_smiles=normalized_smiles,
            input_format=resolved_format,
            unsupported_reason=(
                "The molecule is peptide-like, but it falls outside PepLink v1's reliable reverse scope "
                "(for example sidechain crosslinks, terminal modifications, or non-canonical residues)."
            ),
        )

    return PeptideParseResult(
        sequence=sequence,
        is_cyclic=is_cyclic,
        cyclization="head_to_tail" if is_cyclic else "linear",
        normalized_smiles=normalized_smiles,
        input_format=resolved_format,
        unsupported_reason=None,
    )
