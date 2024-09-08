from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from ase import Atoms
from ase.optimize import BFGS, FIRE, LBFGS, BFGSLineSearch
from sevenn_runner import SevenNetCalculator
import time
import json

max_steps = 50
import gc


def last_non_hydrogen_idx(atomic_nums: np.ndarray):
    return np.argmax(atomic_nums != 1)

# NOTE: this will work poortly for hydrogenated molecules
def calculate_displacement_vector(conf, idx_of_last_atom: int):
    first_atom = conf.GetAtomPosition(0)
    last_atom = conf.GetAtomPosition(idx_of_last_atom)
    return np.array([last_atom.x - first_atom.x, 
                     last_atom.y - first_atom.y, 
                     last_atom.z - first_atom.z])

# uses Rodrigues' rotation formula
def rotation_matrix_to_z(v):
    """
    Compute the rotation matrix that rotates vector v towards the z-axis.
    
    :param v: numpy array, the vector to be rotated
    :return: numpy array, 3x3 rotation matrix
    """
    # Ensure v is a unit vector
    v = v / np.linalg.norm(v)
    
    # Define z-axis
    z = np.array([0, 0, 1])
    
    # Compute rotation axis
    a = np.cross(v, z)
    
    # If v is already aligned with z, return identity matrix
    if np.allclose(a, 0):
        return np.eye(3)
    
    # Normalize rotation axis
    a = a / np.linalg.norm(a)
    
    # Compute rotation angle
    cos_theta = np.dot(v, z)
    cos_theta = cos_theta / 2 # divide by 2, so it's not exactly aligned with z-axis, but in the right direction
    sin_theta = np.sqrt(1 - cos_theta**2)
    
    # Compute rotation matrix using Rodrigues' formula
    K = np.array([
        [0, -a[2], a[1]],
        [a[2], 0, -a[0]],
        [-a[1], a[0], 0]
    ])
    
    R = np.eye(3) + sin_theta * K + (1 - cos_theta) * np.dot(K, K)
    
    return R

def translate_molecule(mol, x=0, y=0, z=0):
    # Make sure the molecule has a conformer
    if mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol)
    
    conf = mol.GetConformer()
    
    # Create translation vector
    translation = np.array([x, y, z])
    
    # Apply translation to each atom
    for atom_idx in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(atom_idx)
        new_pos = np.array([pos.x, pos.y, pos.z]) + translation
        conf.SetAtomPosition(atom_idx, new_pos)
    
    return mol

def rotate_coordinates_3d(coords, rotation_matrix):
    """
    Rotate an array of 3D coordinates using a given 3x3 rotation matrix.
    
    :param coords: numpy array of shape (n, 3) where n is the number of points
    :param rotation_matrix: 3x3 numpy array representing the rotation matrix
    :return: rotated coordinates as a numpy array of shape (n, 3)
    """
    return np.dot(coords, rotation_matrix.T)

def get_molecules(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    idx_of_last_atom = mol.GetNumAtoms() - 1
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)

    # make the molecules face the positive z direction
    direction = calculate_displacement_vector(mol.GetConformer(), idx_of_last_atom)
    r = rotation_matrix_to_z(direction)

    first_atom = mol.GetConformer().GetAtomPosition(0)
    mol = translate_molecule(mol, x=-first_atom.x, y=-first_atom.y, z=-first_atom.z)


    atomic_nums = []
    coords = []
    for i in range(mol.GetNumAtoms()):
        atomic_nums.append(mol.GetAtomWithIdx(i).GetAtomicNum())
        pos = mol.GetConformer().GetAtomPosition(i)
        coords.append([pos.x, pos.y, pos.z])

    coords = np.array(coords)
    atomic_nums = np.array(atomic_nums)
    rotated_coords = rotate_coordinates_3d(coords, r)

    return atomic_nums, rotated_coords, idx_of_last_atom


def to_ase_atoms(atomic_nums: np.ndarray, coords: np.ndarray):
    # lattice_matrix = np.array([
    #     [1000, 0, 0],
    #     [0, 1000, 0],
    #     [0, 0, 1000]
    # ])
    atoms = Atoms(
        numbers=atomic_nums,
        positions=coords,
        # cell=lattice_matrix,
        # pbc=(True, True, True),
        pbc=(False, False, False),
    )
    return atoms

class LoggingBFGS(BFGS):
    def __init__(self, atoms, logfile=None, trajectory=None, coords_log=[]):
        super().__init__(atoms, logfile, trajectory)
        self.coords_log = coords_log
        self.steps_taken = 0

    def log(self):
        super().log()
        print("step: ", len(self.coords_log))
        # Log the fractional positions at each step
        self.coords_log.append(
            self.atoms.get_positions().copy().tolist()
        )

    def step(self, f=None):
        if self.steps_taken >= max_steps:
            return
        # if self.steps_taken % 10 == 0:
        #     gc.collect()
            # Clear Python's internal freelists
            # gc.collect()
            # gc.collect()

            # Clear Python's module cache
            # sys.modules.clear()

        super().step(f)
        self.steps_taken += 1

class LoggingFIRE(FIRE):
    def __init__(self, atoms, logfile=None, trajectory=None, coords_log=[], max_steps=None):
        super().__init__(atoms, logfile=logfile, trajectory=trajectory)
        self.coords_log = coords_log
        self.steps_taken = 0
        self.max_steps = max_steps

    def log(self):
        super().log()
        print("step: ", len(self.coords_log))
        # Log the fractional positions at each step
        self.coords_log.append(
            self.atoms.get_positions().copy().tolist()
        )

    def step(self, f=None):
        if self.max_steps is not None and self.steps_taken >= self.max_steps:
            return
        super().step(f)
        self.steps_taken += 1


def relax(atomic_nums: np.ndarray, coords: np.ndarray):
    sevennet_0_cal = SevenNetCalculator("7net-0", device="auto")  # 7net-0, SevenNet-0, 7net-0_22May2024, 7net-0_11July2024 ...

    # properties = ["energy", "forces", "stress"]

    # set initial positions
    system = to_ase_atoms(atomic_nums=atomic_nums, coords=coords)

    # create the calculator
    print(sevennet_0_cal.device)
    system.calc = sevennet_0_cal

    # dyn = BFGS(system)
    coords_log = []
    start = time.time()
    dyn = LoggingBFGS(system, coords_log=coords_log)
    # dyn = BFGS(system)
    # dyn = LoggingFIRE(system, coords_log=coords_log)
    # dyn = LBFGS(system) 
    # dyn = FIRE(system)
    while not dyn.steps_taken >= max_steps:
        dyn.run(steps=1, fmax=0.0001)
    end = time.time()
    print("relax time: ", end - start)

    return coords_log

def grow_two_molecules(smiles: str):
    mol1_atomic_nums, mol1_coords, mol1_last_non_hydrogen_idx = get_molecules(smiles)
    mol2_atomic_nums, mol2_coords, mol2_last_non_hydrogen_idx = get_molecules(smiles)

    # mol1_non_hydrogen_idx = last_non_hydrogen_idx(mol1_atomic_nums)
    mol1_non_hydrogen_idx = mol1_last_non_hydrogen_idx
    mol2_non_hydrogen_idx = 0

    # move mol2 right next to the last non-hydrogen atom of mol1
    amount_to_move = mol1_coords[mol1_non_hydrogen_idx] - mol2_coords[mol2_non_hydrogen_idx]

    buffer = np.array([3, 3, 3])
    # translate the second molecule to the furthest point of the first molecule
    for i in range(len(mol2_coords)):
        mol2_coords = mol2_coords + amount_to_move + buffer

    print(mol1_coords[-1])
    print(mol2_coords[mol2_last_non_hydrogen_idx+1])
    # remove the last hydrogen of the first molecule
    mol1_atomic_nums = mol1_atomic_nums[:-1]
    mol1_coords = mol1_coords[:-1]
    # remove the first hydrogen of the first molecule
    mol2_atomic_nums = np.delete(mol2_atomic_nums, mol2_last_non_hydrogen_idx +1, axis=0)
    mol2_coords = np.delete(mol2_coords, mol2_last_non_hydrogen_idx+1, axis=0)

    # now we have two molecules near each other. run md to get them to attach

    atomic_nums = np.concatenate((mol1_atomic_nums, mol2_atomic_nums))
    coords = np.concatenate((mol1_coords, mol2_coords))

    coords_log = relax(atomic_nums, coords)

    relaxation = {
        "atomic_nums": atomic_nums.tolist(),
        "coords": coords_log,
    }
    json.dump(relaxation, open("relaxation.json", "w"))