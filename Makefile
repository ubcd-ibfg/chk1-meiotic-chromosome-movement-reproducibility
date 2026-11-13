# Reproducibility package – Camacho-Cabañas et al.
# Usage:  make help
#
# Scripts must run in dependency order (some read what earlier ones write), so the
# stages below are chained explicitly. Do not use `make -j`.

SHELL := bash
.DEFAULT_GOAL := help

ROOT    := $(CURDIR)
CODE    := $(ROOT)/code
VENV    := $(ROOT)/.venv
PYTHON  ?= $(VENV)/bin/python
NPERM   ?= 200          # label permutations; NPERM=20 for a quick check (permutation P values will differ)
export NPERM

# Script 02 variant: "02_series1_classifiers" (default) or "02_series1_classifiers_fixed"
SCRIPT02 ?= 02_series1_classifiers

RUN = cd $(CODE) && $(PYTHON)

.PHONY: help setup setup-pip all all-fixed verify \
        pombe sga cerevisiae suppfig4 \
        s01 s02 s03 s06 s05 s07 s09 s10 s08 s11 s12 s17 s13 s14 s15 s16 \
        clean distclean

help:
	@echo "Setup"
	@echo "  make setup          create .venv with uv (uv sync --locked)"
	@echo "  make setup-pip      create .venv with pip from requirements.txt"
	@echo ""
	@echo "Run"
	@echo "  make all            run every script in dependency order, then verify"
	@echo "  make all-fixed      same, using 02_series1_classifiers_fixed.py"
	@echo "  make pombe          S. pombe imaging   (01 02 03 06 05 07 09 10 08)"
	@echo "  make sga            SGA screen         (11 12 17)"
	@echo "  make cerevisiae     S. cerevisiae      (13 14 15)"
	@echo "  make suppfig4       Supp. Fig. 4       (16)"
	@echo "  make sNN            a single script, e.g. make s07 (needs its inputs)"
	@echo "  make verify         compare the 142 numbers quoted in the paper"
	@echo ""
	@echo "Clean"
	@echo "  make clean          remove __pycache__"
	@echo "  make distclean      also remove output/ (Source_Data/ is kept)"
	@echo ""
	@echo "Variables: NPERM=$(NPERM)  PYTHON=$(PYTHON)  SCRIPT02=$(SCRIPT02)"

# ---------------------------------------------------------------- setup
setup:
	cd $(ROOT) && uv sync --locked

setup-pip:
	python3.11 -m venv $(VENV)
	$(VENV)/bin/pip install -r $(ROOT)/requirements.txt

# ---------------------------------------------------------------- full pipeline
all: pombe sga cerevisiae suppfig4 verify

all-fixed:
	$(MAKE) all SCRIPT02=02_series1_classifiers_fixed

# ---------------------------------------------------------------- stages (chained: order matters)
pombe: s01 s02 s03 s06 s05 s07 s09 s10 s08
sga: s11 s12 s17
cerevisiae: s13 s14 s15
suppfig4: s16

# The explicit prerequisites below enforce the order 01→02→03→06→05→07→09→10→08→11→12→17→13→14→15→16
s01:
	$(RUN) 01_series1_descriptors.py
s02: s01
	$(RUN) $(SCRIPT02).py
s03: s02
	$(RUN) 03_motif_statistics.py
s06: s03
	$(RUN) 06_suppfig2b_ungrouped.py
s05: s06            # 06 must run before 05
	$(RUN) 05_vector_panels.py
s07: s05
	$(RUN) 07_motif_panels.py
s09: s07
	$(RUN) 09_fig1f_vector.py
s10: s09
	$(RUN) 10_suppfig2_composite.py
s08: s10
	$(RUN) 08_source_data.py
s11: s08
	$(RUN) 11_sga_scores.py
s12: s11
	$(RUN) 12_sga_panels.py
s17: s12
	$(RUN) 17_sga_source_data.py
s13: s17
	$(RUN) 13_cerevisiae_analysis.py
s14: s13
	$(RUN) 14_cerevisiae_montage.py
s15: s14
	$(RUN) 15_cerevisiae_source_data.py
s16: s15
	$(RUN) 16_suppfig4_reproducibility.py

# ---------------------------------------------------------------- verification
verify:
	cd $(ROOT) && $(PYTHON) code/verify_paper_numbers.py

# ---------------------------------------------------------------- clean
clean:
	find $(CODE) -name __pycache__ -type d -prune -exec rm -rf {} +

distclean: clean
	rm -rf $(ROOT)/output
