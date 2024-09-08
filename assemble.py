import get_molecules
import json
from sevenn_runner import SevenNetCalculator
import math
import numpy as np


def create_single_chain():
    sevennet_0_cal = SevenNetCalculator("7net-0", device="auto")  # 7net-0, SevenNet-0, 7net-0_22May2024, 7net-0_11July2024 ...
    print(f"running on device {sevennet_0_cal.device}")
    smiles="NCC(=O)NCCCCCC(=O)"
    atomic_nums, coords_log, last_non_hydrogen_idx_on_main_chain = get_molecules.grow_two_molecules(sevennet_0_cal, smiles)

    relax_batches = [{
        "atomic_nums": atomic_nums.tolist(),
        "relax_len": len(coords_log),
    }]

    num_monomers = 5
    for i in range(num_monomers - 2): # -2 since we already have the first 2 monomers
        atomic_nums, coords_log, last_non_hydrogen_idx_on_main_chain = get_molecules.grow_on_chain(sevennet_0_cal, relax_batches[-1]["atomic_nums"], coords_log, last_non_hydrogen_idx_on_main_chain, smiles)
        relax_batches.append({
            "atomic_nums": atomic_nums.tolist(),
            "relax_len": len(coords_log),
        })

    coords_log = coords_log + get_molecules.relax(sevennet_0_cal, relax_batches[-1]["atomic_nums"], coords_log[-1], max_steps=100)
    relax_batches.append({
        "atomic_nums": relax_batches[-1]["atomic_nums"],
        "relax_len": len(coords_log),
    })

    relaxation = {
        "frames": []
    }

    curr_idx = 0
    for relax_batch in relax_batches:
        atomic_nums = relax_batch["atomic_nums"]
        relax_len = relax_batch["relax_len"]
        while curr_idx < relax_len:
            coords = coords_log[curr_idx]
            relaxation["frames"].append({
                "atomic_nums": atomic_nums,
                "coords": coords.tolist(),
            })
            curr_idx += 1
    json.dump(relaxation, open("relaxation.json", "w"), separators=(',', ':'))

def get_initial_coords(num_polymers:int):
    cube_side_len = math.sqrt(num_polymers)
    assert cube_side_len %1 == 0, "Number of num_polymers must be a perfect square"
    initial_positions = []

    distance_between_chains = 10
    cube_side_len = int(cube_side_len)
    for i in range(cube_side_len):
        for j in range(cube_side_len):
            pos = [i*distance_between_chains, j*distance_between_chains, 0]
            initial_positions.append(pos)
    return np.array(initial_positions)

def create_bulk_polymer():
    sevennet_0_cal = SevenNetCalculator("7net-0", device="auto")  # 7net-0, SevenNet-0, 7net-0_22May2024, 7net-0_11July2024 ...
    print(f"running on device {sevennet_0_cal.device}")
    smiles="NCC(=O)NCCCCCC(=O)"
    coords_log = None
    atomic_nums = None
    relax_batches =[]
    for initial_coord in get_initial_coords(4):
        atomic_nums, coords_log, last_non_hydrogen_idx_on_main_chain = get_molecules.grow_two_molecules(sevennet_0_cal, smiles, initial_coord, atomic_nums, coords_log)

        relax_batches.append({
            "atomic_nums": atomic_nums.tolist(),
            "relax_len": len(coords_log),
        })

        num_monomers = 5
        for i in range(num_monomers - 2): # -2 since we already have the first 2 monomers
            atomic_nums, coords_log, last_non_hydrogen_idx_on_main_chain = get_molecules.grow_on_chain(sevennet_0_cal, relax_batches[-1]["atomic_nums"], coords_log, last_non_hydrogen_idx_on_main_chain, smiles)
            relax_batches.append({
                "atomic_nums": atomic_nums.tolist(),
                "relax_len": len(coords_log),
            })

    # final relaxation
    print("running final relaxation")
    coords_log = coords_log + get_molecules.relax(sevennet_0_cal, relax_batches[-1]["atomic_nums"], coords_log[-1], max_steps=50)
    relax_batches.append({
        "atomic_nums": relax_batches[-1]["atomic_nums"],
        "relax_len": len(coords_log),
    })

    relaxation = {
        "frames": []
    }

    curr_idx = 0
    for relax_batch in relax_batches:
        atomic_nums = relax_batch["atomic_nums"]
        relax_len = relax_batch["relax_len"]
        while curr_idx < relax_len:
            coords = coords_log[curr_idx]
            relaxation["frames"].append({
                "atomic_nums": atomic_nums,
                "coords": coords.tolist(),
            })
            curr_idx += 1
    json.dump(relaxation, open("relaxation.json", "w"), separators=(',', ':'))

if __name__ == "__main__":
    create_bulk_polymer()