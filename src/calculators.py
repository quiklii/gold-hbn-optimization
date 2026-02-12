# src/calculators.py

import torch
from ase.calculators.lj import LennardJones
from mace.calculators import mace_mp
from huggingface_hub import hf_hub_download


def detect_device():
    if torch.cuda.is_available():
        return "cuda"
    return "cpu"

def get_calculator(mode="mace", device=None):
    if mode == "mace":
        if device is None:
            device = detect_device()
            print(f"Using device: {device}")
        path = hf_hub_download(repo_id="mace-foundations/mace-mh-1", filename="mace-mh-1.model")
        return mace_mp(model=path, device=device, default_dtype="float64", head="oc20_usemppbe")
    elif mode == 'lj':
        return LennardJones(epsilon=0.02, sigma=2.8)
