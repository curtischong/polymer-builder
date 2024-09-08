import get_molecules
import json

smiles="NCC(=O)NCCCCCC(=O)"
atomic_nums, coords_log = get_molecules.grow_two_molecules(smiles)

relaxation = {
    "atomic_nums": atomic_nums.tolist(),
    "coords": coords_log.tolist(),
}
json.dump(relaxation, open("relaxation.json", "w"), separators=(',', ':'))