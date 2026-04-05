from __future__ import annotations

import json
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from rdkit import RDLogger

from PepLink import (
    UnsupportedPeptideError,
    aa_seqs_to_smiles,
    clear_registered_noncanonical_aas,
    from_dbaasp_record,
    list_supported_noncanonical_aas,
    register_noncanonical_aa,
    register_noncanonical_aas,
    register_noncanonical_aas_from_csv,
    smiles_to_aa_seqs,
)


ROOT = Path(__file__).resolve().parents[1]
DATASET = json.loads((ROOT / "all_peptides_data.json").read_text())
BY_ID = {item["id"]: item for item in DATASET}
RDLogger.DisableLog("rdApp.*")


SUPPORTED_FORWARD_IDS = [11, 151, 157, 10360, 57, 105, 1079, 4419, 6917, 19104, 11913, 14660, 17263, 17273, 19191]


def cyclic_rotations(sequence: str) -> set[str]:
    return {sequence[index:] + sequence[:index] for index in range(len(sequence))}


class AaSeqsToSmilesApiTests(unittest.TestCase):
    def tearDown(self) -> None:
        clear_registered_noncanonical_aas()

    def test_examples_from_dataset(self) -> None:
        for record_id in SUPPORTED_FORWARD_IDS:
            with self.subTest(record_id=record_id):
                kwargs = from_dbaasp_record(BY_ID[record_id]).to_api_kwargs()
                smiles = aa_seqs_to_smiles(**kwargs)
                self.assertIsInstance(smiles, str)
                self.assertTrue(smiles)

    def test_unsupported_bond_type_is_rejected(self) -> None:
        with self.assertRaises(UnsupportedPeptideError):
            aa_seqs_to_smiles(
                "AC",
                intrachain_bonds=[{"position1": 1, "position2": 2, "type": "ETH", "chain_participating": "SSB"}],
            )

    def test_unmapped_unusual_residue_is_rejected(self) -> None:
        with self.assertRaises(UnsupportedPeptideError):
            kwargs = from_dbaasp_record(BY_ID[9851]).to_api_kwargs()
            aa_seqs_to_smiles(**kwargs)

    def test_known_est_edge_case_is_currently_rejected(self) -> None:
        with self.assertRaises(UnsupportedPeptideError):
            aa_seqs_to_smiles(**from_dbaasp_record(BY_ID[5779]).to_api_kwargs())

    def test_registered_noncanonical_aa_is_used_by_default_forward_api(self) -> None:
        register_noncanonical_aa("MyAA", "N[C@@H](CC)C(=O)O")
        smiles = aa_seqs_to_smiles(
            "AXA",
            unusual_amino_acids=[{"position": 2, "name": "MyAA"}],
        )
        self.assertIsInstance(smiles, str)
        self.assertTrue(smiles)

    def test_bulk_registered_noncanonical_aas_are_listed(self) -> None:
        register_noncanonical_aas(
            {
                "MyAA": "N[C@@H](CC)C(=O)O",
                "MyAA2": "N[C@@H](CO)C(=O)O",
            }
        )
        supported = list_supported_noncanonical_aas()
        self.assertIn("MyAA", supported)
        self.assertIn("MyAA2", supported)

    def test_registered_noncanonical_aas_from_csv_are_used_by_forward_api(self) -> None:
        with tempfile.TemporaryDirectory() as target:
            csv_path = Path(target) / "custom_noncanonical.csv"
            csv_path.write_text("name,SMILES\nMyAA,N[C@@H](CC)C(=O)O\n", encoding="utf-8")
            register_noncanonical_aas_from_csv(csv_path)
        smiles = aa_seqs_to_smiles(
            "AXA",
            unusual_amino_acids=[{"position": 2, "name": "MyAA"}],
        )
        self.assertIsInstance(smiles, str)
        self.assertTrue(smiles)

    def test_multimer_and_coordination_records_are_rejected(self) -> None:
        with self.assertRaises(UnsupportedPeptideError):
            from_dbaasp_record(BY_ID[1])
        with self.assertRaises(UnsupportedPeptideError):
            from_dbaasp_record(BY_ID[15])


class SmilesToAaSeqsApiTests(unittest.TestCase):
    def test_linear_round_trip(self) -> None:
        smiles = aa_seqs_to_smiles(**from_dbaasp_record(BY_ID[11]).to_api_kwargs())
        parsed = smiles_to_aa_seqs(smiles)
        self.assertEqual(parsed.sequence, BY_ID[11]["sequence"])
        self.assertFalse(parsed.is_cyclic)
        self.assertEqual(parsed.cyclization, "linear")
        self.assertIsNone(parsed.unsupported_reason)

    def test_head_to_tail_round_trip(self) -> None:
        smiles = aa_seqs_to_smiles(**from_dbaasp_record(BY_ID[105]).to_api_kwargs())
        parsed = smiles_to_aa_seqs(smiles)
        self.assertIn(parsed.sequence, cyclic_rotations(BY_ID[105]["sequence"]))
        self.assertTrue(parsed.is_cyclic)
        self.assertEqual(parsed.cyclization, "head_to_tail")
        self.assertIsNone(parsed.unsupported_reason)

    def test_complex_aa_seqs_to_smiles_case_is_marked_unsupported_in_reverse(self) -> None:
        smiles = aa_seqs_to_smiles(**from_dbaasp_record(BY_ID[57]).to_api_kwargs())
        parsed = smiles_to_aa_seqs(smiles)
        self.assertIsNone(parsed.sequence)
        self.assertIsNotNone(parsed.unsupported_reason)

    def test_invalid_smiles_returns_structured_failure(self) -> None:
        parsed = smiles_to_aa_seqs("not-a-smiles")
        self.assertIsNone(parsed.sequence)
        self.assertIsNotNone(parsed.unsupported_reason)

    @unittest.skipUnless(importlib.util.find_spec("selfies") is not None, "selfies is not installed in the local test environment")
    def test_selfies_round_trip(self) -> None:
        selfies_text = aa_seqs_to_smiles("AC", output_format="selfies")
        parsed = smiles_to_aa_seqs(selfies_text, input_format="selfies")
        self.assertEqual(parsed.sequence, "AC")
        self.assertIsNone(parsed.unsupported_reason)


class PackagingSmokeTests(unittest.TestCase):
    def test_local_install_smoke(self) -> None:
        with tempfile.TemporaryDirectory() as target:
            subprocess.run(
                [sys.executable, "-m", "pip", "install", ".", "--no-deps", "--no-build-isolation", "-t", target],
                cwd=ROOT,
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            code = (
                "import sys;"
                f"sys.path.insert(0, {target!r});"
                "from PepLink import aa_seqs_to_smiles;"
                "print(bool(aa_seqs_to_smiles('AC')))"
            )
            completed = subprocess.run(
                [sys.executable, "-c", code],
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.stdout.strip(), "True")
