# read the relaxation.json file and get the last frame
import json
import numpy as np

from get_molecules import to_ase_atoms
from sevenn_batch import SevennBatchEval

with open("relaxation.json", "r") as f:
    relaxation = json.load(f)

last_frame = relaxation["frames"][-1]

atomic_nums = last_frame["atomic_nums"]
coords = last_frame["coords"]

num_stretch_frames = 20
total_stretch_amount = 0.1 # 10%

all_atoms = []

for i in range(num_stretch_frames):
    current_stretch_amount = total_stretch_amount * i / num_stretch_frames
    stretch_amount = np.array([1, 1, 1 + current_stretch_amount])
    coords = stretch_amount * coords
    system = to_ase_atoms(atomic_nums=atomic_nums, coords=coords)
    all_atoms.append(system)

sevenn_batch_eval = SevennBatchEval(
    batch_size=1  # yes. we need to set it to this low or the gpu will run out of memory
)

results = sevenn_batch_eval.eval(all_atoms)

frames = []
for i in range(len(results)):
    res = results[i]
    coords = all_atoms[i].get_positions().tolist()
    energy = np.clip(res["inferred_total_energy"], a_min=-9999, a_max=9999).tolist()
    forces = np.clip(res["inferred_force"], a_min=-9999, a_max=9999).tolist()
    stress = np.clip(res["inferred_stress"], a_min=-9999, a_max=9999).tolist()

    frames.append({
        "atomic_nums": atomic_nums,
        "coords": coords,
        "energy": energy,
        "forces": forces,
        "stress": stress,
    })

pull = {
    "frames": frames
}
json.dump(pull, open("pull.json", "w"), separators=(',', ':'))