# src/geometry.py

from typing import Tuple, Optional, List
import numpy as np
from ase import Atoms
from ase.build import make_supercell
from ase.constraints import FixAtoms
from ase.geometry import get_distances


def create_hbn_slab(size: Tuple[int, int, int] = (6, 6, 1),
                    vacuum: float = 20.0,
                    constraint_mode: str = 'all') -> Atoms:
    a = 2.51
    c = vacuum

    cell = [
        [a, 0.0, 0.0],
        [a/2, a * np.sqrt(3)/2, 0.0],
        [0.0, 0.0, c]
    ]

    slab = Atoms(symbols=['B', 'N'], cell=cell, pbc=[True, True, False])
    slab.set_scaled_positions([[0.0, 0.0, 0.5], [1/3, 2/3, 0.5]])

    P = [[size[0], 0, 0],
         [0, size[1], 0],
         [0, 0, size[2]]]
    slab = make_supercell(slab, P)

    if constraint_mode == 'all':
        constraint = FixAtoms(indices=[atom.index for atom in slab])
        slab.set_constraint(constraint)
    elif constraint_mode == 'anchor':
        constraint = FixAtoms(indices=[0])
        slab.set_constraint(constraint)

    return slab


def get_z_range(slab: Atoms, height_range: Tuple[float, float] = (2.0, 5.0)) -> Tuple[float, float]:
    positions = slab.get_positions()
    z_surface = positions[:, 2].max()
    return (z_surface + height_range[0], z_surface + height_range[1])


def add_random_cluster(slab: Atoms,
                       n_atoms: int,
                       element: str = 'Au',
                       seed: Optional[int] = None,
                       min_distance: float = 2.0,
                       max_attempts: int = 1000) -> Atoms:
    rng = np.random.default_rng(seed)

    system = slab.copy()

    z_min, z_max = get_z_range(slab)

    cluster_atoms = Atoms()

    cell = system.get_cell().array

    margin_frac = 0.05

    for i in range(n_atoms):
        placed = False

        for attempt in range(max_attempts):
            u = rng.uniform(margin_frac, 1.0 - margin_frac)
            v = rng.uniform(margin_frac, 1.0 - margin_frac)

            z = rng.uniform(z_min, z_max)

            pos_xy = u * cell[0] + v * cell[1]
            pos_candidate = np.array([pos_xy[0], pos_xy[1], z])
            _, dists_slab = get_distances(
                pos_candidate, system.positions, cell=cell, pbc=system.pbc)

            if len(cluster_atoms) > 0:
                _, dists_cluster = get_distances(
                    pos_candidate, cluster_atoms.positions, cell=cell, pbc=system.pbc)
                min_d = min(np.min(dists_slab), np.min(dists_cluster))
            else:
                min_d = np.min(dists_slab)

            if min_d >= min_distance:
                atom = Atoms(element, positions=[pos_candidate])
                cluster_atoms.extend(atom)
                placed = True
                break

        if not placed:
            raise RuntimeError(
                f"Nie udało się umieścić atomu {i+1}/{n_atoms} po {max_attempts} próbach. "
                f"Zwiększ obszar (box) lub zmniejsz min_distance."
            )

    final_system = system + cluster_atoms
    final_system.set_cell(slab.get_cell())
    final_system.set_pbc(slab.get_pbc())

    if slab.constraints:
        final_system.set_constraint(slab.constraints)

    return final_system
