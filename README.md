<img src="https://raw.githubusercontent.com/quiklii/gold-hbn-optimization/refs/heads/master/public/readme_icon.png" width="50" />

# 🧬 GA-Au-hBN: Optymalizacja Klastra Au na h-BN

Repozytorium robocze w ramach laboratorium naukowego.

**Cel:** Znalezienie najstabilniejszej geometrii klastra $Au_6$ na powierzchni heksagonalnego azotku boru (h-BN) przy użyciu Algorytmu Genetycznego (GA) w środowisku ASE.

## 📂 Struktura Projektu

```text
gold-hbn-optimization/
├── data/               # Pliki wynikowe: bazy danych (.db), trajektorie (.traj)
├── notebooks/          # Jupyter Notebooks do ewentualnej analizy wykresów i wizualizacji
├── src/                # Kod źródłowy modułów
│   ├── __init__.py
│   ├── calculators.py  # Logika wyboru kalkulatora
│   ├── geometry.py     # Definicja superkomórki h-BN, constraints, boxa
│   └── ga_setup.py     # Konfiguracja operatorów genetycznych
├── run.py              # Główny punkt wejścia (CLI)
├── pyproject.toml      # Zależności projektu
└── README.md           # Dokumentacja
```
