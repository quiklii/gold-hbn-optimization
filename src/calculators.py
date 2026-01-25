# src/calculators.py

from ase.calculators.lj import LennardJones

def get_calculator(mode="lj", **kwargs):
    mode = (mode or "lj").lower().strip()
    if mode == "lj":
        return LennardJones(
            epsilon=kwargs.get("epsilon", 0.02),  # eV
            sigma=kwargs.get("sigma", 2.8),      # A
            rc=kwargs.get("rc", 8.0),            # A
            smooth=kwargs.get("smooth", True)
        )
    raise ValueError(f"Nieznany tryb: {mode}")
