from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
from rdkit.Chem import Atom


def get_molecules(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)

    atomic_nums = []
    coords = []
    for i in range(mol.GetNumAtoms()):
        atomic_nums.append(mol.GetAtomWithIdx(i).GetAtomicNum())
        pos = mol.GetConformer().GetAtomPosition(i)
        coords.append([pos.x, pos.y, pos.z])

    return atomic_nums, coords



def view_molecule(mol, width=800, height=600, style=None, surface=None):
    """
    Generate an interactive 3D viewer for an RDKit molecule.
    """
    pdb = Chem.MolToPDBBlock(mol)
    viewer = py3Dmol.view(width=width, height=height)
    viewer.addModel(pdb, "pdb")
    if style is None:
        style = {"stick": {"radius": 0.15}, "sphere": {"radius": 0.3}}
    viewer.setStyle(style)
    if surface:
        viewer.addSurface(py3Dmol.VDW, surface)
    viewer.zoomTo()