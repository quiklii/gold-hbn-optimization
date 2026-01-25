# src/geometry.py

"""
    Moduł do tworzenia geometrii
"""

from ase.atoms import Atoms
from ase.build import make_supercell
from ase.constraints import FixAtoms
import numpy as np

def create_hbn_slab(size=(6, 6, 1), vacuum=20.0):
    a = 2.504
    c=vacuum

    # komorka heksagonalna
    cell = [
        [a, 0.0, 0.0],
        [a/2, a * np.sqrt(3)/2, 0.0],
        [0.0, 0.0, c]
    ]

    # ustawienia atomow
    slab = Atoms(symbols=['B', 'N'], cell=cell, pbc=[True, True, False])
    slab.set_scaled_positions([
        [0.0, 0.0, 1/2], # B
        [1/3, 2/3, 1/2]  # N
    ])

    # tworzenie superkomorki
    P = [[size[0], 0, 0],
         [0, size[1], 0],
         [0, 0, size[2]]]
    slab = make_supercell(slab, P)

    # mrozenie podloza
    constraint = FixAtoms(indices=[atom.index for atom in slab])
    slab.set_constraint(constraint)

    return slab

def get_cluster_box(slab, height_range=(2.0, 5.0), margin=1.0):
    cell = slab.get_cell()
    positions = slab.get_positions()

    z_surface = positions[:, 2].max()

    x_max = cell[0, 0] - margin
    y_max = cell[1, 1] - margin

    box = {
        'x' : (margin, x_max),
        'y' : (margin, y_max),
        'z' : (z_surface + height_range[0], z_surface + height_range[1])
    }

    return box

def add_random_cluster(slab, n_atoms, element='Au', box=None, seed=None):
    rng = np.random.default_rng(seed)

    if box is None:
        box = get_cluster_box(slab)

    # losowe pozycje w pudelku
    positions = np.array([
        rng.uniform(box['x'][0], box['x'][1], n_atoms),
        rng.uniform(box['y'][0], box['y'][1], n_atoms),
        rng.uniform(box['z'][0], box['z'][1], n_atoms)
    ]).T

    # klaster
    cluster = Atoms([element] * n_atoms, positions=positions)

    #laczymy z podlazem
    system = slab + cluster
    system.set_cell(slab.get_cell())
    system.set_pbc(slab.get_pbc())

    system.set_constraint(FixAtoms(indices=list(range(len(slab)))))

    return system


# === Testy ===
if __name__ == "__main__":
    print("🧪 Test 1: Tworzenie h-BN...")
    slab = create_hbn_slab(size=(6, 6, 1), vacuum=20.0)
    print(f"   ✅ {len(slab)} atomów")
    print(f"   Stosunek B:N = {slab.get_chemical_symbols().count('B')}:{slab.get_chemical_symbols().count('N')}")
    print(f"   Komórka XY: {slab.get_cell()[0, 0]:.2f} x {slab.get_cell()[1, 1]:.2f} Å")

    print("\n🧪 Test 2: Pudełko dla Au...")
    box = get_cluster_box(slab, height_range=(2.0, 5.0))
    print(f"   ✅ X: {box['x']}")
    print(f"   ✅ Y: {box['y']}")
    print(f"   ✅ Z: {box['z']}")

    print("\n🧪 Test 3: Dodanie Au₆...")
    system = add_random_cluster(slab, n_atoms=6, element='Au', seed=42)
    print(f"   ✅ Całkowicie: {len(system)} atomów (h-BN: {len(slab)}, Au: {len(system) - len(slab)})")
    print(f"   Au pozycje Z: {system.get_positions()[len(slab):, 2]}")