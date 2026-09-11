"""Guards for the residue definitions corrected in 0.1.3.

Twenty-six shipped structures were wrong: some named a different compound
entirely, some were written as free zwitterions that cannot form a peptide
bond, and some were multi-component salts. These tests pin the corrections so a
future edit cannot quietly restore them.
"""

from __future__ import annotations

import unittest

from rdkit import Chem, RDLogger
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from PepLink import aa_seqs_to_smiles
from PepLink.data import load_default_aa_smiles

RDLogger.DisableLog("rdApp.*")

# code -> formula the DBAASP record states for that residue.
CORRECTED_FORMULAS = {
    "Allo-Thr": "C4H9NO3",
    "ARGol": "C6H16N4O",
    "CCA": "C7H12O2",
    "I3CA": "C9H7NO2",
    "Iaa": "C5H6N2O2",
    "Iac": "C6H6N2O2",
    "LYS-(CH3)3": "C9H21N2O2+",
    "Nspe": "C10H13NO2",
}

# Residues that must be usable inside a chain, not only at the N-terminus.
INTERNAL_RESIDUES = ("Allo-Thr", "LYS-(CH3)3", "MIM", "Aic", "Agb", "Nae", "3-Me-Trp")


class CorrectedDefinitions(unittest.TestCase):
    def setUp(self) -> None:
        self.definitions = load_default_aa_smiles()

    def test_every_definition_parses(self) -> None:
        for code, smiles in self.definitions.items():
            with self.subTest(code=code):
                self.assertIsNotNone(
                    Chem.MolFromSmiles(smiles), f"{code} does not parse: {smiles}"
                )

    def test_corrected_formulas_match_the_source_records(self) -> None:
        for code, expected in CORRECTED_FORMULAS.items():
            with self.subTest(code=code):
                molecule = Chem.MolFromSmiles(self.definitions[code])
                self.assertEqual(CalcMolFormula(molecule), expected)

    def test_only_the_known_suspect_remains_a_multi_component_salt(self) -> None:
        """A residue is one connected structure; a salt cannot be polymerised.

        ``LAP`` stays as shipped because its source record contradicts itself:
        the stated formula C16H34N2O4 matches neither reading of the name it
        gives, so replacing it would be a guess. Any other salt is a regression.
        """

        offenders = sorted(
            code for code, smiles in self.definitions.items() if "." in smiles
        )
        self.assertEqual(offenders, ["LAP"])

    def test_internal_residues_build_inside_a_chain(self) -> None:
        """A zwitterionic or non-polymerisable definition fails this."""

        for code in INTERNAL_RESIDUES:
            with self.subTest(code=code):
                smiles = aa_seqs_to_smiles(
                    "AXA",
                    unusual_amino_acids=[
                        {"position": 2, "modificationType": {"name": code}}
                    ],
                )
                self.assertIsNotNone(Chem.MolFromSmiles(smiles))

    def test_iaa_and_iac_are_no_longer_the_same_compound(self) -> None:
        """Both shipped indole-3-acetic acid; they are imidazole compounds."""

        self.assertNotEqual(self.definitions["Iaa"], self.definitions["Iac"])
        for code in ("Iaa", "Iac"):
            self.assertNotIn("c2ccccc2", Chem.MolToSmiles(
                Chem.MolFromSmiles(self.definitions[code])
            ))


if __name__ == "__main__":
    unittest.main()
