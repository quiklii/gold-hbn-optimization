# 🧬 GA-Au-hBN: Optymalizacja Klastra Au na h-BN

Repozytorium robocze w ramach laboratorium naukowego.
**Cel:** Znalezienie najstabilniejszej geometrii klastra $Au_6$ na powierzchni heksagonalnego azotku boru (h-BN) przy użyciu Algorytmu Genetycznego (GA) w środowisku ASE.

---

## 📅 Roadmapa (Plan Rozwoju)

### Faza 1: Setup Środowiska i Fizyka
- [ ] **Struktura CLI:** Skrypt `run.py` obsługujący argumenty (`argparse`).
- [ ] **Geometria (Slab):**
    - [ ] Zdefiniować superkomórkę h-BN ($6\times6$).
    - [ ] Dodać próżnię (~20 Å) w osi Z.
    - [ ] Nałożyć `FixAtoms` na atomy podłoża.
- [ ] **Kalkulatory (Logika w `src/calculators.py`):**
    - [ ] **Hybryda (Default):** Implementacja `SumCalculator` (EMT dla Au + LJ dla Au-BN). *Priorytet.*
    - [ ] **MACE (Feature):** Implementacja obsługi modelu MLIP na GPU.

### Faza 2: Algorytm Genetyczny (Implementacja)
- [ ] **Populacja Startowa:**
    - [ ] Generator 20 losowych klastrów Au₆ w pudełku 2-5 Å nad h-BN.
    - [ ] Weryfikacja `blmin` (2.5 Å dla Au-Au).
- [ ] **Komparator:**
    - [ ] `InteratomicDistanceComparator` (odległość < 0.5 Å = duplikat).
- [ ] **Selekcja:**
    - [ ] Tournament selection (rozmiar turnieju: 3).
- [ ] **Operatory:**
    - [ ] `CutAndSplicePairing` (prawdopodobieństwo: 0.5).
    - [ ] `RattleMutation` (stdev=0.8, prawdopodobieństwo: 0.3).
    - [ ] `PermutationMutation` (jeśli brak MirrorMutation).
    - [ ] `RotationalMutation` (sprawdzić nazwę!).
- [ ] **Pętla Główna:**
    - [ ] Relaksacja BFGS (fmax=0.05 eV/Å, max_steps=200).
    - [ ] Zapis do `gadb.db` (co 1 krok).
    - [ ] Logging (energia best/worst, czas kroku).

### Faza 3: Wykonanie i Raport (przykładowy plan)
- [ ] **Test lokalny:** Puszczenie 10 kroków na CPU (dla sprawdzenia błędów).
- [ ] **Produkcja:** Puszczenie pełnej pętli na RTX (MACE) lub CPU (Hybryda).
- [ ] **Analiza:**
    - [ ] Wyciągnięcie Top 3 struktur.
    - [ ] Wykres zbieżności (Energia vs Krok).
    - [ ] Ocena struktury (2D vs 3D) (?)

---

## 📂 Struktura Projektu

```text
gold-hbn-optimization/
├── data/               # Pliki wynikowe: bazy danych (.db), trajektorie (.traj)
│                       # (Folder ignorowany przez git)
├── notebooks/          # Jupyter Notebooks do ewentualnej analizy wykresów i wizualizacji
├── src/                # Kod źródłowy modułów
│   ├── __init__.py
│   ├── calculators.py  # Logika wyboru kalkulatora (MACE vs EMT+LJ)
│   ├── geometry.py     # Definicja superkomórki h-BN, constraints, boxa
│   └── ga_setup.py     # Konfiguracja operatorów genetycznych
├── run.py              # Główny punkt wejścia (CLI)
├── pyproject.toml      # Zależności projektu
└── README.md           # Dokumentacja i plan