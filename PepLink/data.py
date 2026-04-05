from __future__ import annotations

import csv
from pathlib import Path
from functools import lru_cache
from importlib.resources import files
from threading import RLock
from typing import Dict, Mapping

from ._chem import normalize_smiles
from .exceptions import ValidationError
from .models import STANDARD_AA_CODES_WITH_D


_CANONICAL_AA_NAMES = STANDARD_AA_CODES_WITH_D | {"g"}
_CUSTOM_NONCANONICAL_AA_SMILES: Dict[str, str] = {}
_CUSTOM_NONCANONICAL_AA_LOCK = RLock()


def _load_csv_mapping(relative_path: str) -> Dict[str, str]:
    resource = files("PepLink").joinpath(relative_path)
    with resource.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return {row["name" if "name" in row else "aa"]: row["SMILES"] for row in reader}


@lru_cache(maxsize=1)
def load_default_aa_smiles() -> Dict[str, str]:
    mapping = _load_csv_mapping("data/all_aa_smiles_new_handcrafted.csv")
    mapping.setdefault("g", mapping["G"])
    # DAB is mis-entered in the legacy CSV; normalize it to L-2,4-diaminobutyric acid.
    mapping["DAB"] = "N[C@@H](CCN)C(=O)O"
    return mapping


def _normalize_noncanonical_name(name: str) -> str:
    normalized = name.strip()
    if not normalized:
        raise ValidationError("non-canonical amino-acid name must not be empty")
    if normalized in _CANONICAL_AA_NAMES:
        raise ValidationError(f"{normalized!r} is a reserved canonical amino-acid code and cannot be re-registered")
    return normalized


def _normalize_noncanonical_smiles(name: str, smiles: str) -> str:
    normalized_input = smiles.strip()
    if not normalized_input:
        raise ValidationError(f"SMILES for non-canonical amino acid {name!r} must not be empty")
    try:
        return normalize_smiles(normalized_input)
    except ValueError as exc:
        raise ValidationError(f"Invalid SMILES for non-canonical amino acid {name!r}: {smiles}") from exc


def _validated_noncanonical_mapping(mapping: Mapping[str, str]) -> Dict[str, str]:
    validated: Dict[str, str] = {}
    for name, smiles in mapping.items():
        normalized_name = _normalize_noncanonical_name(name)
        validated[normalized_name] = _normalize_noncanonical_smiles(normalized_name, smiles)
    return validated


def register_noncanonical_aa(name: str, smiles: str) -> str:
    normalized_name = _normalize_noncanonical_name(name)
    normalized_smiles = _normalize_noncanonical_smiles(normalized_name, smiles)
    with _CUSTOM_NONCANONICAL_AA_LOCK:
        _CUSTOM_NONCANONICAL_AA_SMILES[normalized_name] = normalized_smiles
    return normalized_smiles


def register_noncanonical_aas(mapping: Mapping[str, str]) -> Dict[str, str]:
    registered = _validated_noncanonical_mapping(mapping)
    with _CUSTOM_NONCANONICAL_AA_LOCK:
        _CUSTOM_NONCANONICAL_AA_SMILES.update(registered)
    return registered


def clear_registered_noncanonical_aas() -> None:
    with _CUSTOM_NONCANONICAL_AA_LOCK:
        _CUSTOM_NONCANONICAL_AA_SMILES.clear()


def load_noncanonical_aas_from_csv(csv_path: str | Path) -> Dict[str, str]:
    path = Path(csv_path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        name_field = "name" if "name" in fieldnames else "aa" if "aa" in fieldnames else None
        if name_field is None or "SMILES" not in fieldnames:
            raise ValidationError("non-canonical amino-acid CSV must contain columns 'name' (or 'aa') and 'SMILES'")
        raw_mapping: Dict[str, str] = {}
        for line_number, row in enumerate(reader, start=2):
            raw_name = row.get(name_field, "")
            raw_smiles = row.get("SMILES", "")
            if not raw_name.strip():
                raise ValidationError(f"non-canonical amino-acid CSV row {line_number} is missing a residue name")
            if not raw_smiles.strip():
                raise ValidationError(
                    f"non-canonical amino-acid CSV row {line_number} is missing a SMILES value for {raw_name!r}"
                )
            raw_mapping[raw_name] = raw_smiles
    return _validated_noncanonical_mapping(raw_mapping)


def register_noncanonical_aas_from_csv(csv_path: str | Path) -> Dict[str, str]:
    registered = load_noncanonical_aas_from_csv(csv_path)
    with _CUSTOM_NONCANONICAL_AA_LOCK:
        _CUSTOM_NONCANONICAL_AA_SMILES.update(registered)
    return registered


def load_registered_aa_smiles() -> Dict[str, str]:
    mapping = dict(load_default_aa_smiles())
    with _CUSTOM_NONCANONICAL_AA_LOCK:
        mapping.update(_CUSTOM_NONCANONICAL_AA_SMILES)
    return mapping


def list_supported_noncanonical_aas(*, include_custom: bool = True) -> Dict[str, str]:
    if include_custom:
        mapping = load_registered_aa_smiles()
    else:
        mapping = load_default_aa_smiles()
    return {name: mapping[name] for name in sorted(mapping) if name not in _CANONICAL_AA_NAMES}


@lru_cache(maxsize=1)
def load_default_n_terminal_smiles() -> Dict[str, str]:
    return _load_csv_mapping("data/terminal_modifications/n_terminal_smiles_from_PubChem_handcrafted.csv")


@lru_cache(maxsize=1)
def load_default_c_terminal_smiles() -> Dict[str, str]:
    return _load_csv_mapping("data/terminal_modifications/c_terminal_smiles_from_PubChem_handcrafted.csv")


def merged_mapping(defaults: Mapping[str, str], overrides: Mapping[str, str] | None) -> Dict[str, str]:
    merged = dict(defaults)
    if overrides:
        merged.update(overrides)
    return merged
