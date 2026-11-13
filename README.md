# chk1-meiotic-chromosome-movement-reproducibility
### Camacho-Cabañas et al. – *DNA damage signalling reprogrammes meiotic chromosome movement through Chk1*

This folder contains everything needed to regenerate every quantitative panel, the Source Data and the Series 1 rows of Supplementary Table 1 from the original data. The raw data are unmodified.

## Folder layout

| Folder | Contents |
|---|---|
| `raw_data/S_pombe_imaging/` | `horsetail_data.csv`: frame-level ChroMo output, one row per segmented frame of the chromosome mass. `submission_data_AFA_2026/`: per-strain motif detections, the official motif counts per strain (`*_motifs_and_windows_summary.csv`, J. Pinto-Cruz) and the strain-code correspondence. `SuppFig1a_cell_size.xlsx`: meiocyte size from transmitted-light images. |
| `raw_data/SGA_screen/` | `data_heatmap_completo_replica123_sergio_log2fold.csv`: original SGA colony-size table (3,403 deletions × medium × query × replicate). `annotation/`: PomBase GAF (2026-05-21) and GO (go-basic, 2026-07-26) used for the gene categories. |
| `raw_data/S_cerevisiae/` | Original TEL4L-GFP trajectories (46 control and 137 *rdh54Δ* cells) and the two representative movies (Fig. 5F). |
| `code/` | All analysis scripts, `run_all.sh` (runs them in dependency order; table below), `run_all_fixed.sh` (same, with the fixed script 02) and `verify_paper_numbers.py` (checks the outputs against the paper). |
| `original_scripts/` | The lab's original scripts, unmodified: R scripts of the SGA screen (reference only) and `S_cerevisiae/cerevisiae_ML_analysis.py`, whose feature definitions script 13 imports, so this folder **is required**. |
| `Makefile` | Entry point: `make setup`, `make all`, `make verify`, per-dataset targets (see "How to run"). |
| `pyproject.toml`, `uv.lock`, `.python-version` | uv project files: Python 3.11.7 and every package pinned exactly (`uv sync --locked`). `uv.lock` is cross-platform (Linux, macOS, Windows). |
| `requirements.txt` | The same exact versions for pip, including indirect dependencies. |
| `output/` | Tables and figure panels written by the code, one subfolder per dataset. |
| `Source_Data/` | Source Data workbooks, one per dataset; each has a README sheet. |

## Installation

### Requirements

| Requirement | Notes |
|---|---|
| Python **3.11** | Checked with 3.11.7. `pyproject.toml` declares `requires-python = "==3.11.*"`, so uv refuses other versions. |
| [uv](https://docs.astral.sh/uv/) | Recommended installer (uses `uv.lock`). Not needed if you use pip. |
| GNU `make` | Optional; only for the `make` commands below. `bash run_all.sh` works without it. |

Exact package versions (identical in `uv.lock` and `requirements.txt`): numpy 2.4.6, pandas 3.1.1, scipy 1.17.1, scikit-learn 1.9.1, matplotlib 3.11.2, networkx 3.6.1, joblib 1.6.0, openpyxl 3.1.5 and tifffile 2026.3.3.

**scikit-learn 1.9.1 is required for Supp. Fig. 2.** The cross-validation folds of `StratifiedGroupKFold(shuffle=True)` depend on the scikit-learn version. With 1.9.1, script 02 reproduces the published values exactly. With 1.7.2 or older, it stops with an error (see below), and the fixed script gives different P values. Other versions can also change the Fig. 5B/5C permutation P values and feature ranking.

### Option A: uv (recommended)

```bash
cd Reproducibility
make setup                 # = uv sync --locked: creates .venv with Python 3.11.7 and the exact versions of uv.lock
# without make:  uv sync --locked
```

### Option B: pip

```bash
cd Reproducibility
make setup-pip             # = python3.11 -m venv .venv && .venv/bin/pip install -r requirements.txt
# without make:  python3.11 -m venv .venv && .venv/bin/pip install -r requirements.txt
```

Both options create `.venv/` next to `code/`. No activation is needed: the Makefile and the `run_all*.sh` runners call `.venv/bin/python` directly. To use another interpreter, set `PYTHON=/path/to/python`.

## How to run

### With make (recommended)

Run all commands from the `Reproducibility/` folder (where the `Makefile` is).

```bash
make               # same as make help: lists all targets
make all           # every script in dependency order, then verify (~10 min on 88 cores; much longer on a laptop, mostly script 02)
make all NPERM=20  # quick check with 20 permutations (permutation P values will differ)
make verify        # compare the 142 values quoted in the paper with the outputs
```

| Target | What it does |
|---|---|
| `make setup` / `make setup-pip` | Create `.venv` (see Installation). |
| `make all` | Runs the 16 scripts in order, then `verify`. |
| `make all-fixed` | Same, with `02_series1_classifiers_fixed.py` instead of script 02. |
| `make pombe` | *S. pombe* imaging: 01 → 02 → 03 → 06 → 05 → 07 → 09 → 10 → 08. |
| `make sga` | SGA screen: 11 → 12 → 17. Needs the outputs of `make pombe`. |
| `make cerevisiae` | *S. cerevisiae*: 13 → 14 → 15. |
| `make suppfig4` | Supp. Fig. 4: 16. Needs the outputs of the previous stages. |
| `make s01` … `make s17` | One script (e.g. `make s07`). It **also runs every earlier script in the chain**, because each target depends on the previous one. To run a single script without its predecessors, use the Python command in the next section. |
| `make verify` | Runs `code/verify_paper_numbers.py`. |
| `make clean` | Removes `code/__pycache__`. |
| `make distclean` | Also removes `output/`. `Source_Data/` is kept; it is regenerated (overwritten) by the next full run. |

Variables, set on the command line (`make all NPERM=20`):

| Variable | Default | Meaning |
|---|---|---|
| `NPERM` | `200` | Number of label permutations in the permutation tests. |
| `PYTHON` | `.venv/bin/python` | Interpreter used for every script. |
| `SCRIPT02` | `02_series1_classifiers` | Script 02 variant; `make all-fixed` sets it to `02_series1_classifiers_fixed`. |

Do not use `make -j`: the stages read each other's outputs and must run one after another. Stages are not skipped when their outputs exist; every `make` call re-runs what it names.

### With the shell scripts

```bash
cd code
bash run_all.sh            # every script in dependency order
NPERM=20 bash run_all.sh   # quick check with 20 permutations
bash run_all_fixed.sh      # same, with 02_series1_classifiers_fixed.py
cd .. && .venv/bin/python code/verify_paper_numbers.py
```

`run_all_fixed.sh` is the same as `run_all.sh` but uses `02_series1_classifiers_fixed.py`. With scikit-learn 1.9.1 the two give identical results; the fixed script only matters with older versions.

`run_all.sh` and `run_all_fixed.sh` run every script with `../.venv/bin/python`. Set `PYTHON=/path/to/python` to use another interpreter. If `.venv` is missing they fall back to `python3` and print a warning (the Makefile does not fall back: it fails if the interpreter is missing).

### One script on its own

```bash
cd code && ../.venv/bin/python 07_motif_panels.py
```

The scripts create `output/` and `Source_Data/` if they do not exist. The scripts must run in this order, because some read what earlier ones write:
01 → 02 → 03 → 06 → 05 → 07 → 09 → 10 → 08 → 11 → 12 → 17 → 13 → 14 → 15 → 16.
In particular, 06 must run before 05. Each script can also be run on its own once its inputs exist, from any directory.

Folder layout note: `code/` also contains `run_all.sh`/`run_all_fixed.sh`; the `Makefile` sits at the top of the package, next to `pyproject.toml`.

## Scripts and the panels they produce

| Script | Produces | Figure panels |
|---|---|---|
| `common.py` | Shared settings and paths (imported by 01–10) | — |
| `01_series1_descriptors.py` | Per-cell descriptors, statistics vs control, z-scores, PCA, meiocyte size, Supp. Table 1 | Fig. 1C, E, F, G; Supp. Fig. 1A–H; Supp. Table 1 |
| `02_series1_classifiers.py` | Classifier benchmark, random forest vs control (grouped by cell, 200 permutations), Gini importance, 91 pairwise comparisons. Requires scikit-learn 1.9.1: with ≤ 1.7.2 it stops with `ValueError: Found array with 0 sample(s)` (see below) | Supp. Fig. 2A–D |
| `02_series1_classifiers_fixed.py` | Same as 02, but skips cross-validation folds in which train or test has only one class (none occur with scikit-learn 1.9.1, so results are identical there) | Supp. Fig. 2A–D |
| `03_motif_statistics.py` | Motif A/C rates vs the control decomposed in the same batch, exact conditional test, exact 95 % CI, BH-adjusted P, observation windows | numbers for Fig. 2E, J–L; 3J, K; 4D–G |
| `05_vector_panels.py` | Vector panels | Fig. 1G; Supp. Fig. 2A–D |
| `06_suppfig2b_ungrouped.py` | Cross-validation not grouped by cell | Supp. Fig. 2B (grey ticks) |
| `07_motif_panels.py` | Motif panels (Fig. 3K includes the *rdh54Δ hrs1-T16A* vs *rdh54Δ* bracket) | Fig. 2E, 2J–L; Fig. 3J, K; Fig. 4D, E (latrunculin A); Fig. 4F, G (*myo51Δ*) |
| `08_source_data.py` | `Source_Data/Source_Data_S_pombe_imaging.xlsx` | — |
| `09_fig1f_vector.py` | z-score heat map, vector | Fig. 1F |
| `10_suppfig2_composite.py` | Supp. Fig. 2 composite | Supp. Fig. 2 |
| `11_sga_scores.py` | SGA interaction scores (rdh54Δ − control query), ranks, GO-category tests; reproduces the original top-200 list (199/200) | Fig. 4A–C, ranks in the text |
| `12_sga_panels.py` | SGA panels | Fig. 4A, B, C |
| `13_cerevisiae_analysis.py` | *S. cerevisiae* features, balanced-accuracy classifiers with permutation tests, feature importance, ACF, large-displacement episodes | Fig. 5A–E, G, H |
| `14_cerevisiae_montage.py` | Montage of the representative cells | Fig. 5F |
| `15_cerevisiae_source_data.py` | `Source_Data/Source_Data_Fig5_cerevisiae.xlsx` | — |
| `17_sga_source_data.py` | `Source_Data/Source_Data_Fig4_SGA.xlsx` | — |
| `16_suppfig4_reproducibility.py` | Motif C in *rdh54Δ* across the four series, pooled rate ratio, observation windows (windows per cell, events per window); `Source_Data/Source_Data_SuppFig4.xlsx` | Supp. Fig. 4A, B |
| `verify_paper_numbers.py` | Compares 142 numbers quoted in the paper with the outputs, within rounding (run from `Reproducibility/`) | — |

`make all` and `run_all.sh` run everything in the right order. Note: `output/S_pombe_imaging/suppfig3_pairwise.csv` holds Supp. Fig. 2D (the file keeps an earlier name).


## Statistical conventions

- **Unit of analysis:** the meiocyte. Descriptors are averaged per cell before testing; motif frequency is expressed as events per cell.
- **Descriptors:** two-sided Mann–Whitney U test vs control; Benjamini–Hochberg across the 13 genotypes.
- **Classifiers (*S. pombe*):** five-fold cross-validation grouped by cell, with balanced folds; significance from 200 label permutations at the cell level.
- **Motifs:** each strain is compared with the control decomposed in the same batch, using an exact conditional test for two Poisson rates with exact 95 % CI.
- **SGA:** score = mean(*rdh54Δ* replicates) − mean(control replicates), separately on YES and YES + 15 mM HU. Missing colonies are kept as missing, not set to 0.
- **S. cerevisiae:** pixel size 0.0663 µm, 1 frame per minute, balanced accuracy over repeated stratified 5-fold cross-validation, 200 permutations.

## Experimental series (*S. pombe*)

| Series | Imaging system | Summary file (`raw_data/S_pombe_imaging/submission_data_AFA_2026/`) | Figures |
|---|---|---|---|
| 1 | Olympus IX81 / Yokogawa spinning disk | `Figure_1/` and `Supplementary_Figure_1/` | Fig. 1, 2E; Supp. Figs 1–3 |
| 2 | DeltaVision | `Figure_rdh54_double_mutants/`, `Figure_Microtubules_Actin/` (microtubules) | Fig. 2J–L, 3J, K |
| 3 | Nikon Ti2 E / Andor Dragonfly | `Figure_Rodrigo/` | Fig. 4F, G |
| 4 | Nikon Ti2 E / Andor Dragonfly | `Figure_Microtubules_Actin/` (latrunculin A) | Fig. 4D, E |



## Citation

If you use this code or data, please cite it as below (see also [`CITATION.cff`](CITATION.cff)):

> Pinto Cruz, J., & Fernández-Álvarez, A. (2026). *chk1-meiotic-chromosome-movement-reproducibility* (v0.1.1) [Computer software]. https://github.com/ubcd-ibfg/chk1-meiotic-chromosome-movement-reproducibility

```bibtex
@software{chk1_reproducibility_2026,
  author  = {Pinto Cruz, Jesús and Fernández-Álvarez, Alfonso},
  title   = {chk1-meiotic-chromosome-movement-reproducibility},
  version = {0.1.1},
  year    = {2026},
  url     = {https://github.com/ubcd-ibfg/chk1-meiotic-chromosome-movement-reproducibility}
}
```
