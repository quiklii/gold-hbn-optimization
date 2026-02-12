import argparse
from pathlib import Path
from datetime import datetime

import numpy as np
from ase.optimize import BFGS
from ase.io import write
from ase.constraints import FixAtoms

from src.geometry import create_hbn_slab, add_random_cluster
from src.calculators import get_calculator
from src.ga_setup import configure_ga


def set_raw_score(atoms):
    energy = atoms.get_potential_energy()
    if 'key_value_pairs' not in atoms.info:
        atoms.info['key_value_pairs'] = {}
    atoms.info['key_value_pairs']['raw_score'] = - \
        energy

    if 'data' not in atoms.info:
        atoms.info['data'] = {}


def parse_args():
    p = argparse.ArgumentParser(description="GA: Au6 on h-BN")
    p.add_argument("--db-name", default="gadb.db")
    p.add_argument("--n-initial", type=int, default=20)
    p.add_argument("--n-steps", type=int, default=100)
    p.add_argument("--population-size", type=int, default=20)
    p.add_argument("--fmax", type=float, default=0.05)
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--mode", type=str, default="mace", choices=["mace", "lj"])
    return p.parse_args()


def relax(atoms, calc, fmax=0.05, steps=200):
    atoms.calc = calc
    try:
        BFGS(atoms, logfile=None).run(fmax=fmax, steps=steps)
        forces = atoms.get_forces()
        converged = np.max(np.linalg.norm(forces, axis=1)) < fmax
        return converged, atoms.get_potential_energy()
    except Exception as e:
        print(f"  Relax error: {e}")
        return False, None


def main():
    args = parse_args()

    slab = create_hbn_slab(size=(6, 6, 1), vacuum=20.0)
    slab_len = len(slab)
    calc = get_calculator(mode=args.mode)

    slab.calc = calc
    E_slab = slab.get_potential_energy()
    print(
        f"Slab energy: {E_slab:.4f} eV for {len(slab)} atoms")
    # ----------------------------------------------------------

    ga = configure_ga(slab, args.db_name, n_Au=6,
                      population_size=args.population_size)
    population, selector, db = ga["population"], ga["selector"], ga["db"]

    population.update()
    if len(population.pop) == 0:
        for i in range(args.n_initial):
            system = add_random_cluster(
                slab, n_atoms=6, element="Au", seed=args.seed + i)
            ok, e_total = relax(system, calc, args.fmax)

            if ok:
                e_relative = e_total - E_slab
                set_raw_score(system)
                system.info['key_value_pairs']['raw_score'] = -e_relative

                db.add_relaxed_candidate(system)
                print(
                    f"  [{i + 1}/{args.n_initial}] E_ads = {e_relative:.4f} eV")
            else:
                print(f"  [{i + 1}/{args.n_initial}] FAIL")

    for step in range(args.n_steps):
        population.update()

        if len(population.pop) < 2:
            print(f"[{step + 1}] Not enough candidates, skip")
            continue

        parents = population.get_two_candidates()
        offspring, desc = selector.get_new_individual(parents)

        if offspring is None:
            print(f"[{step + 1}] Offspring=None, skip")
            continue

        offspring.set_cell(slab.get_cell())
        offspring.set_pbc(slab.get_pbc())
        offspring.set_constraint(FixAtoms(indices=list(range(slab_len))))

        ok, e_total = relax(offspring, calc, args.fmax)

        if ok:
            e_relative = e_total - E_slab

            set_raw_score(offspring)
            offspring.info['key_value_pairs']['raw_score'] = -e_relative

            db.add_relaxed_candidate(offspring)

            best_total = min(a.get_potential_energy() for a in population.pop)
            best_rel = best_total - E_slab

            print(
                f"[{step + 1}] {desc}: E_ads={e_relative:.4f} | best={best_rel:.4f}")
        else:
            print(f"[{step + 1}] {desc}: FAIL")

    population.update()
    if population.pop:
        best = min(population.pop, key=lambda a: a.get_potential_energy())

        final_e = best.get_potential_energy() - E_slab

        out = Path("data") / \
            f"best_{datetime.now():%Y%m%d_%H%M%S}_E{final_e:.2f}.xyz"
        write(out, best)
        print(f"\nBest found: E_ads = {final_e:.4f} eV -> {out}")


if __name__ == "__main__":
    main()
