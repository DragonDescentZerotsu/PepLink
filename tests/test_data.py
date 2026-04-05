from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from rdkit import RDLogger

from PepLink.data import (
    clear_registered_noncanonical_aas,
    list_supported_noncanonical_aas,
    load_noncanonical_aas_from_csv,
    load_default_aa_smiles,
    load_default_c_terminal_smiles,
    load_default_n_terminal_smiles,
    register_noncanonical_aa,
    register_noncanonical_aas_from_csv,
)


ROOT = Path(__file__).resolve().parents[1]
RDLogger.DisableLog("rdApp.*")


class DataCoverageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.dataset = json.loads((ROOT / "all_peptides_data.json").read_text())

    def tearDown(self) -> None:
        clear_registered_noncanonical_aas()

    def test_noncanonical_mapping_count(self) -> None:
        mapping = load_default_aa_smiles()
        canonical = set("ACDEFGHIKLMNPQRSTVWY") | {aa.lower() for aa in "ACDEFGHIKLMNPQRSTVWY"} | {"g"}
        noncanonical = {name for name in mapping if name not in canonical}
        self.assertEqual(len(noncanonical), 420)

    def test_public_noncanonical_listing_excludes_canonical_codes(self) -> None:
        mapping = list_supported_noncanonical_aas(include_custom=False)
        self.assertEqual(len(mapping), 420)
        self.assertNotIn("A", mapping)
        self.assertNotIn("a", mapping)
        self.assertIn("1-NAL", mapping)

    def test_registered_noncanonical_smiles_is_normalized_and_listed(self) -> None:
        normalized = register_noncanonical_aa("MyAA", "OC(=O)[C@H](N)CC")
        mapping = list_supported_noncanonical_aas()
        self.assertEqual(mapping["MyAA"], normalized)

    def test_load_noncanonical_aas_from_csv_returns_normalized_mapping(self) -> None:
        with tempfile.TemporaryDirectory() as target:
            csv_path = Path(target) / "custom_noncanonical.csv"
            csv_path.write_text("name,SMILES\nMyAA,OC(=O)[C@H](N)CC\n", encoding="utf-8")
            mapping = load_noncanonical_aas_from_csv(csv_path)
        self.assertIn("MyAA", mapping)
        self.assertEqual(mapping["MyAA"], "CC[C@@H](N)C(=O)O")

    def test_register_noncanonical_aas_from_csv_updates_runtime_registry(self) -> None:
        with tempfile.TemporaryDirectory() as target:
            csv_path = Path(target) / "custom_noncanonical.csv"
            csv_path.write_text("name,SMILES\nMyAA,N[C@@H](CC)C(=O)O\n", encoding="utf-8")
            register_noncanonical_aas_from_csv(csv_path)
        mapping = list_supported_noncanonical_aas()
        self.assertIn("MyAA", mapping)

    def test_terminal_modification_coverage_matches_dataset(self) -> None:
        n_mapping = load_default_n_terminal_smiles()
        c_mapping = load_default_c_terminal_smiles()
        n_terms = {((item.get("nTerminus") or {}).get("name")) for item in self.dataset if item.get("nTerminus")}
        c_terms = {((item.get("cTerminus") or {}).get("name")) for item in self.dataset if item.get("cTerminus")}
        self.assertEqual(n_terms, set(n_mapping))
        self.assertEqual(c_terms, set(c_mapping))
