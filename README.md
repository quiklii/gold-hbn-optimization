<p style="text-align:center;">
  <img src="https://raw.githubusercontent.com/quiklii/gold-hbn-optimization/refs/heads/master/public/readme.png"  />
</p>

# 🧬 GA-Au-hBN: Au Cluster Optimization on h-BN using MACE potential

Working repository for a scientific lab project.

**Goal:** Find the most stable geometry of an $Au_6$ cluster adsorbed on a hexagonal boron nitride (h-BN) surface using a **Genetic Algorithm (GA)** in the **ASE** environment.

---

## ▶️ How to run

### 1) Create and activate a virtual environment
```bash
python -m venv .venv

# Windows (PowerShell)
.venv\Scripts\Activate.ps1

# Linux/macOS
source .venv/bin/activate
```
### 2) Install dependencies (from pyproject.toml)
Install the project in editable mode:
```bash
pip install -e .
```
For CUDA, do not rely on install PyTorch first from the official selector [here](https://pytorch.org/get-started/locally/).

### 3) Run the GA optimization
```bash
python run.py
```
For more options, check:
```bash
python run.py --help
```
Outputs (database + best structures) are written to the `data/` directory.

---

## 🧠 How it works

This repository solves a **global structure optimization** problem: finding the lowest-energy geometry of an **Au₆** cluster adsorbed on an **h-BN (α-BN)** surface.  
The search is done with a **steady-state Genetic Algorithm (GA)** in **ASE**, while energies/forces are evaluated using the **MACE** machine-learning interatomic potential.

### System model
- **Substrate:** h-BN slab in a **6×6 supercell** with **periodic boundary conditions**.
- **Fixed surface:** all B and N atoms are constrained using `FixAtoms` (rigid substrate approximation).
- **Adsorbate:** only the **Au₆** atoms are optimized (free to move above the surface).

### Energy model (MACE)
- Calculator: **MACE-MH1**
- Head: **`oc20_usemppbe`** (recommended for surface–adsorbate systems)

### Optimization strategy: steady-state GA (ASE)
Unlike generational GAs (where the whole population is replaced each generation), this project uses a **steady-state** approach: the population is updated continuously by inserting improved individuals and discarding worse ones.

**Workflow**
1. **Initialization**
   - Generate an initial population of physically valid Au₆ geometries above the slab.
   - Locally relax each candidate with **BFGS**.
2. **Evolution loop**
   - Repeat for a fixed number of steps:
     - **Parent selection:** tournament selection (ASE `Population` implementation).
     - **Variation:** create an offspring using genetic operators:
       - crossover / recombination (geometric)
       - mutation (random displacements, reflections, rotations)
     - **Local relaxation:** relax the offspring (BFGS).
     - **Insert/update:** add the new relaxed structure to the database and update the population, maintaining diversity while favoring low-energy candidates.
3. **Output**
   - The run stores energies and best structures in `data/` (DB + best geometries).

### Fitness / objective
The GA minimizes the **adsorption energy**, i.e. it searches for the most thermodynamically stable adsorbed configuration (most negative \(E_{ads}\)).

## 📈 Results

**Best structure example:**
<p style="text-align:center;">
  <img src="https://raw.githubusercontent.com/quiklii/gold-hbn-optimization/refs/heads/master/public/best_struct.png"/>
</p>

**Short description:**

- Final evolutionary optimization used a population of **50** individuals over **2000** steps.
- The full run took **~11 h** and converged to a stable **Au₆** cluster.
- Best configuration reached an adsorption energy of **E_ads = −10.76 eV** (≈ **−1.79 eV/atom**), found early (~**100** steps) and remained stable for the rest of the run (deep minimum on the PES).
- All Au atoms lie within **z = 13.90–14.15 Å** (spread **0.25 Å**), indicating a **planar** cluster geometry (consistent with literature).
- The average distance to the h-BN layer (at **z = 10.00 Å**) is **~4.0 Å**, suggesting adsorption dominated by **van der Waals** interactions.

**Why an (approximately) equilateral-triangle shape?**  
The optimized Au₆ cluster tends to form a highly symmetric, compact 2D arrangement because this **maximizes the number of Au–Au bonds** (higher coordination) while keeping **bond lengths and forces as uniform as possible**. High symmetry often correlates with a **low-energy minimum** on the potential energy surface: fewer distorted bonds → lower strain energy. On weakly interacting substrates like h-BN (vdW-dominated), the substrate imposes only a small perturbation, so the cluster geometry is governed mainly by **internal Au–Au bonding**, favoring a near-equilateral, close-packed motif.