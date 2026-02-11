# src/calculators.py

from ase.calculators.lj import LennardJones
from mace.calculators import mace_mp


def get_calculator(mode="mace", device="cpu"):
    if mode == "mace":
        return mace_mp(model="medium", device=device, default_dtype="float32")
    return LennardJones(epsilon=0.02, sigma=2.8)
