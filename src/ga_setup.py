# src/ga_setup.py

from ase_ga.data import DataConnection, PrepareDB
from ase_ga.population import Population
from ase_ga.standard_comparators import InteratomicDistanceComparator
from ase_ga.utilities import closest_distances_generator
from ase_ga.cutandsplicepairing import CutAndSplicePairing
from ase_ga.standardmutations import MirrorMutation, RattleMutation, RotationalMutation
from ase_ga.offspring_creator import OperationSelector

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent


def configure_ga(slab, filename, n_Au=6, population_size=100):
    blmin = closest_distances_generator(
        atom_numbers=[79, 5, 7],
        ratio_of_covalent_radii=0.8
    )

    data_dir = PROJECT_DIR / 'data'
    data_dir.mkdir(exist_ok=True)

    db_path = data_dir / filename

    if not db_path.exists():
        print(f"Creating new database: {db_path}")
        PrepareDB(
            db_file_name=str(db_path),
            stoichiometry=[79] * n_Au,
            cell=slab.get_cell(),
            pbc=slab.get_pbc()
        )

    db = DataConnection(db_path)

    comparator = InteratomicDistanceComparator(
        n_top=n_Au,
        pair_cor_cum_diff=0.015,
        pair_cor_max=0.7,
        dE=0.02,
        mic=True
    )

    population = Population(
        data_connection=db,
        population_size=population_size,
        comparator=comparator
    )

    crossover = CutAndSplicePairing(
        slab=slab,
        n_top=n_Au,
        blmin=blmin
    )
    rattle = RattleMutation(
        blmin=blmin,
        n_top=n_Au
    )
    mirror = MirrorMutation(
        blmin=blmin,
        n_top=n_Au
    )
    rotation = RotationalMutation(
        blmin=blmin,
        n_top=n_Au
    )

    selector = OperationSelector(
        [0.3, 0.3, 0.2, 0.2],
        [crossover, rattle, mirror, rotation]
    )

    return {
        'population': population,
        'selector': selector,
        'db': db,
        'comparator': comparator
    }
