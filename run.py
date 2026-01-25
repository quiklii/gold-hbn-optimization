# run.py
"""
Runner GA: Au6 na h-BN (ASE + ase-ga) z kalkulatorem LJ i relaksacją BFGS.

Wymaga:
- src/geometry.py: create_hbn_slab, add_random_cluster
- src/calculators.py: get_calculator (LJ)
- src/ga_setup.py: configure_ga (ase_ga + PrepareDB + operatory)
"""

import argparse
import sys
from pathlib import Path
from datetime import datetime

import numpy as np
from ase.optimize import BFGS
from ase.io import write
from ase.constraints import FixAtoms
from ase.calculators.singlepoint import SinglePointCalculator

from src.geometry import create_hbn_slab, add_random_cluster
from src.calculators import get_calculator
from src.ga_setup import configure_ga


def parse_args():
    p = argparse.ArgumentParser(description="GA: Au6 on h-BN (LJ + BFGS, ase-ga)")
    p.add_argument("--db-name", type=str, default="gadb.db", help="Nazwa pliku DB w folderze data/")
    p.add_argument("--n-initial", type=int, default=20, help="Ile struktur startowych (seed)")
    p.add_argument("--n-steps", type=int, default=100, help="Liczba kroków GA")
    p.add_argument("--population-size", type=int, default=20, help="Rozmiar populacji")
    p.add_argument("--fmax", type=float, default=0.05, help="Próg zbieżności BFGS [eV/Å]")
    p.add_argument("--relax-steps", type=int, default=200, help="Maks. liczba kroków BFGS")
    p.add_argument("--seed", type=int, default=7, help="Seed RNG")
    # LJ params (pasują do src/calculators.py)
    p.add_argument("--epsilon", type=float, default=0.02, help="LJ epsilon [eV]")
    p.add_argument("--sigma", type=float, default=2.8, help="LJ sigma [Å]")
    p.add_argument("--rc", type=float, default=8.0, help="LJ cutoff [Å]")
    p.add_argument("--no-smooth", action="store_true", help="Wyłącz smooth w LJ")
    return p.parse_args()


def ensure_slab_frozen(atoms, slab_len: int):
    """Zamraża pierwsze slab_len atomów (podłoże), reszta (Au) jest ruchoma."""
    atoms.set_constraint(FixAtoms(indices=list(range(slab_len))))


def relax_and_attach_singlepoint(atoms, calculator, fmax=0.05, max_steps=200):
    """
    Relaksuje atoms BFGS.
    Sprawdza zbieżność przez max normę sił (wersjo-odporne, nie używa opt.converged()).
    Po relaksacji podpina SinglePointCalculator (energy + forces).
    """
    atoms.calc = calculator

    try:
        opt = BFGS(atoms, logfile=None)
        opt.run(fmax=fmax, steps=max_steps)
    except Exception as e:
        return False, None, f"relax error: {e}"

    try:
        forces = atoms.get_forces()
        fmax_now = float(np.max(np.linalg.norm(forces, axis=1)))
        if fmax_now > fmax:
            return False, None, f"not converged (fmax_now={fmax_now:.3f} > fmax={fmax:.3f})"

        e = float(atoms.get_potential_energy())
    except Exception as e:
        return False, None, f"energy/forces error: {e}"

    atoms.calc = SinglePointCalculator(atoms, energy=e, forces=forces)
    atoms.info["energy"] = e
    return True, e, "ok"


def db_num_relaxed(db):
    """
    Różne wersje ase-ga mają różne nazwy liczników.
    Próbujemy kilku opcji.
    """
    for name in (
        "get_number_of_relaxed_candidates",
        "get_number_of_relaxed_structures",
        "get_number_of_relaxed_steps",
        "get_number_of_relaxed",
    ):
        if hasattr(db, name):
            try:
                return int(getattr(db, name)())
            except Exception:
                pass
    # fallback: jeśli nie ma licznika, spróbujmy pobrać populację później (zwróć None)
    return None


def db_add_relaxed(db, atoms, description: str):
    """
    Adapter do zapisu relaxed kandydata do DB.
    Najczęściej w ase-ga jest add_relaxed_candidate(atoms, ...)
    """
    if hasattr(db, "add_relaxed_candidate"):
        # typowy ase-ga
        return db.add_relaxed_candidate(atoms, description=description)

    if hasattr(db, "add_relaxed_step"):
        # jeśli ktoś ma taki interfejs (starsze/zmodyfikowane)
        return db.add_relaxed_step(atoms, description=description)

    raise AttributeError(
        "Nie znalazłem metody zapisu relaxed w DataConnection. "
        "Sprawdź dir(db) i dopasuj db_add_relaxed()."
    )


def seed_initial_population(slab, db, calculator, n_initial: int, fmax: float, relax_steps: int, seed: int, n_Au: int):
    print(f"\nSeeding: generuję {n_initial} losowych struktur startowych...")
    slab_len = len(slab)
    accepted = 0

    for i in range(n_initial):
        system = add_random_cluster(slab, n_atoms=n_Au, element="Au", seed=seed + i)
        # bezpieczeństwo: zawsze zamroź slab w systemie (operators czasem gubią constraints)
        ensure_slab_frozen(system, slab_len)

        ok, e, msg = relax_and_attach_singlepoint(system, calculator, fmax=fmax, max_steps=relax_steps)
        if not ok:
            print(f"  [{i+1:>2}/{n_initial}] FAIL ({msg})")
            continue

        db_add_relaxed(db, system, description=f"initial_{i}")
        accepted += 1
        print(f"  [{i+1:>2}/{n_initial}] OK  E={e:.6f} eV")

    print(f"Seed done: accepted {accepted}/{n_initial}\n")


def run(args):
    print("=" * 70)
    print("GA: Au6 on h-BN (LJ + BFGS, ase-ga)")
    print("=" * 70)
    print(f"DB: data/{args.db_name}")
    print(f"n_initial={args.n_initial} | n_steps={args.n_steps} | pop_size={args.population_size}")
    print(f"BFGS: fmax={args.fmax} eV/Å | relax_steps={args.relax_steps}")
    print(f"LJ: epsilon={args.epsilon} eV | sigma={args.sigma} Å | rc={args.rc} Å | smooth={not args.no_smooth}")
    print("=" * 70)

    # 1) Slab
    slab = create_hbn_slab(size=(6, 6, 1), vacuum=20.0)
    slab_len = len(slab)
    n_Au = 6  # w tym projekcie na razie stałe Au6

    # 2) Calculator (LJ)
    calculator = get_calculator(
        mode="lj",
    )

    # 3) GA setup (DB + Population + Selector)
    ga = configure_ga(
        slab=slab,
        filename=args.db_name,
        n_Au=n_Au,
        population_size=args.population_size,
    )
    population = ga["population"]
    selector = ga["selector"]
    db = ga["db"]

    # 4) Seed (tylko jeśli DB pusta)
    n_relaxed = db_num_relaxed(db)
    if n_relaxed is None:
        # fallback: spróbujmy zupdate’ować populację i zobaczyć czy coś jest
        try:
            population.update()
            current = getattr(population, "pop", [])
            n_relaxed = len(current)
        except Exception:
            n_relaxed = 0

    if n_relaxed == 0:
        seed_initial_population(
            slab=slab,
            db=db,
            calculator=calculator,
            n_initial=args.n_initial,
            fmax=args.fmax,
            relax_steps=args.relax_steps,
            seed=args.seed,
            n_Au=n_Au,
        )
    else:
        print(f"\nDB już ma relaxed candidates: {n_relaxed} -> pomijam seeding.\n")

    # 5) Main GA loop
    print("Start GA loop...")
    for step in range(args.n_steps):
        print(f"\n[STEP {step+1}/{args.n_steps}]")
        try:
            population.update()
        except Exception as e:
            print(f"population.update() failed: {e}")
            continue

        # rodzice z populacji
        try:
            parents = population.get_two_candidates()
            # zwykle to jest (a1, a2)
            if isinstance(parents, (list, tuple)) and len(parents) == 2:
                a1, a2 = parents[0], parents[1]
            else:
                print("Błąd, skip")
                continue
        except Exception as e:
            print(f"Nie mogę pobrać rodziców: {e}")
            continue

        # wybór operatora + offspring
        try:
            offspring, desc = selector.get_new_individual([a1, a2])
        except Exception as e:
            print(f"Operator failed: {e}")
            continue

        if offspring is None:
            print("offspring=None -> skip")
            continue

        # zabezpieczenia geometrii: cell/pbc/constraints
        offspring.set_cell(slab.get_cell())
        offspring.set_pbc(slab.get_pbc())
        ensure_slab_frozen(offspring, slab_len)

        ok, e, msg = relax_and_attach_singlepoint(offspring, calculator, fmax=args.fmax, max_steps=args.relax_steps)
        if not ok:
            print(f"Relax FAIL ({msg})")
            continue

        print(f"Operator: {desc} | E={e:.6f} eV")
        try:
            db_add_relaxed(db, offspring, description=f"step_{step}_{desc}")
        except Exception as e2:
            print(f"DB write FAIL: {e2}")

        # szybki podgląd populacji (jeśli pop ma Atoms-y z info["energy"])
        try:
            pop_list = getattr(population, "pop", [])
            Es = [a.info.get("energy") for a in pop_list if hasattr(a, "info") and "energy" in a.info]
            if Es:
                print(f"Population: best={min(Es):.6f} eV | worst={max(Es):.6f} eV | n={len(Es)}")
        except Exception:
            pass

    # 6) Save best (z aktualnej populacji)
    try:
        population.update()
    except Exception:
        pass

    pop_list = getattr(population, "pop", [])
    if not pop_list:
        print("\nBrak struktur w populacji na koniec.")
        return

    def energy_of(a):
        if hasattr(a, "info") and "energy" in a.info:
            return float(a.info["energy"])
        try:
            return float(a.get_potential_energy())
        except Exception:
            return float("inf")

    best = min(pop_list, key=energy_of)
    best_E = energy_of(best)

    out_dir = Path("data")
    out_dir.mkdir(exist_ok=True)
    out_file = out_dir / f"best_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xyz"
    write(out_file, best)
    print("\n" + "=" * 70)
    print(f"DONE. Best E = {best_E:.6f} eV")
    print(f"Saved: {out_file}")
    print("=" * 70)


def main():
    args = parse_args()
    try:
        run(args)
    except KeyboardInterrupt:
        print("\nPrzerwano (Ctrl+C)")
        sys.exit(1)
    except Exception as e:
        print(f"\nBłąd krytyczny: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

