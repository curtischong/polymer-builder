import get_molecules
import json

smiles="NCC(=O)NCCCCCC(=O)"
atomic_nums, coords_log, last_non_hydrogen_idx_on_main_chain = get_molecules.grow_two_molecules(smiles)

num_monomers = 5
for i in range(num_monomers):
    atomic_nums, coords_log, last_non_hydrogen_idx_on_main_chain = get_molecules.grow_on_chain(atomic_nums, coords_log, last_non_hydrogen_idx_on_main_chain, smiles)

relaxation = {
    "atomic_nums": atomic_nums.tolist(),
    "coords": [cl.tolist() for cl in coords_log],
}
json.dump(relaxation, open("relaxation.json", "w"), separators=(',', ':'))