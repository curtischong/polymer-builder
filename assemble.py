import get_molecules
import json
from sevenn_runner import SevenNetCalculator


sevennet_0_cal = SevenNetCalculator("7net-0", device="auto")  # 7net-0, SevenNet-0, 7net-0_22May2024, 7net-0_11July2024 ...
smiles="NCC(=O)NCCCCCC(=O)"
atomic_nums, coords_log, last_non_hydrogen_idx_on_main_chain = get_molecules.grow_two_molecules(sevennet_0_cal, smiles)

num_monomers = 5
for i in range(num_monomers - 1):
    atomic_nums, coords_log, last_non_hydrogen_idx_on_main_chain = get_molecules.grow_on_chain(sevennet_0_cal, atomic_nums, coords_log, last_non_hydrogen_idx_on_main_chain, smiles)

relaxation = {
    "atomic_nums": atomic_nums.tolist(),
    "coords": [cl.tolist() for cl in coords_log],
}
json.dump(relaxation, open("relaxation.json", "w"), separators=(',', ':'))