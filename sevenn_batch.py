import logging

import ase
import numpy as np
import sevenn._keys as KEY
import torch
from ase.calculators.singlepoint import SinglePointCalculator
from sevenn._const import LossType
from sevenn.atom_graph_data import AtomGraphData
from sevenn.train.dataload import graph_build
from sevenn.train.dataset import AtomGraphDataset
from sevenn.util import (
    AverageNumber,
    model_from_checkpoint,
    postprocess_output,
    pretrained_name_to_path,
    squared_error,
    to_atom_graph_list,
)
from torch_geometric.loader import DataLoader

logger = logging.getLogger(__name__)


class SevennBatchEval:
    def __init__(self, batch_size: int):
        self.batch_size = batch_size
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        logger.info(f"Using device {self.device}")

        model_name = "7net-0"
        checkpoint = pretrained_name_to_path(model_name)
        self.model, config = model_from_checkpoint(checkpoint)
        self.cutoff = config[KEY.CUTOFF]
        self.type_map = config[KEY.TYPE_MAP]

        self.model.to(self.device)
        self.model.set_is_batch_data(True)
        self.model.eval()

    def eval(self, all_atoms: list[ase.Atoms]) -> list[AtomGraphData]:
        num_cores = 4

        # Prepare atoms list (same as before)
        atoms_list = self._prepare_atoms_list(all_atoms)

        data_list = graph_build(atoms_list, self.cutoff, num_cores=num_cores)
        inference_set = AtomGraphDataset(data_list, self.cutoff)

        inference_set.x_to_one_hot_idx(self.type_map)
        inference_set.toggle_requires_grad_of_data(KEY.POS, True)

        loss_types = [LossType.ENERGY, LossType.FORCE, LossType.STRESS]
        l2_err = {k: AverageNumber() for k in loss_types}
        infer_list = inference_set.to_list()

        try:
            # important that shuffle=False!!! (so the material IDs do not get mixed up)
            loader = DataLoader(infer_list, batch_size=self.batch_size, shuffle=False)
            output_list = []

            for batch in loader:
                # Move entire batch to device and ensure consistent dtype
                batch = batch.to(self.device, non_blocking=True)

                output = self.model(batch)
                output = output.detach().cpu()  # Move output to CPU immediately
                results = postprocess_output(output, loss_types)
                for loss_type in loss_types:
                    l2_err[loss_type].update(squared_error(*results[loss_type]))
                output_list.extend(to_atom_graph_list(output))

                # Clear CUDA cache after each batch
                if self.device.type == "cuda":
                    torch.cuda.empty_cache()

        except Exception as e:
            logger.warn(e)
            logger.info("Keeping 'info' failed. Try with separated info")
            infer_list, _info_list = inference_set.seperate_info()
            loader = DataLoader(infer_list, batch_size=self.batch_size, shuffle=False)
            output_list = []

            for batch in loader:
                # Move entire batch to device and ensure consistent dtype
                batch = batch.to(self.device, non_blocking=True)

                output = self.model(batch)
                output = output.detach().cpu()  # Move output to CPU immediately
                results = postprocess_output(output, loss_types)
                for loss_type in loss_types:
                    l2_err[loss_type].update(squared_error(*results[loss_type]))
                output_list.extend(to_atom_graph_list(output))

                # Clear CUDA cache after each batch
                if self.device.type == "cuda":
                    torch.cuda.empty_cache()

        return output_list

    def _prepare_atoms_list(self, all_atoms):
        # Logic from poscars_to_atoms (same as before)
        stress_dummy = np.array([0, 0, 0, 0, 0, 0])
        calc_results = {"energy": 0, "free_energy": 0, "stress": stress_dummy}
        atoms_list = []
        for atoms in all_atoms:
            natoms = len(atoms.get_atomic_numbers())
            dummy_force = np.zeros((natoms, 3))
            dummy_calc_res = calc_results.copy()
            dummy_calc_res["forces"] = dummy_force
            atoms = SinglePointCalculator(atoms, **dummy_calc_res).get_atoms()
            atoms_list.append(atoms)
        return atoms_list
