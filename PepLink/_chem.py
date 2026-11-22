from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize


if hasattr(Chem.Atom, "GetValence"):
    def _explicit_valence(atom: Chem.Atom) -> int:
        return atom.GetValence(Chem.ValenceType.EXPLICIT)
else:
    # Older RDKit releases expose only the dedicated valence getters.
    _explicit_valence = Chem.Atom.GetExplicitValence


def normalize_mol(mol: Chem.Mol) -> Chem.Mol:
    """Normalize charge/state so comparisons are stable across equivalent inputs."""
    cleaned = rdMolStandardize.Cleanup(Chem.Mol(mol))
    parent = rdMolStandardize.FragmentParent(cleaned)
    uncharger = rdMolStandardize.Uncharger()
    return uncharger.uncharge(parent)


def normalize_smiles(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return canonical_smiles(mol)


def canonical_smiles(mol: Chem.Mol, *, kekule: bool = False) -> str:
    return Chem.MolToSmiles(
        normalize_mol(mol),
        canonical=True,
        isomericSmiles=True,
        kekuleSmiles=kekule,
    )
