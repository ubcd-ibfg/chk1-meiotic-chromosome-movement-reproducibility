#!/usr/bin/env bash
# Regenerates every output and Source Data workbook from raw_data/, in dependency order.
# Usage (from this folder):  bash run_all.sh        (NPERM=200 by default; NPERM=20 bash run_all.sh for a quick check)
set -euo pipefail
cd "$(dirname "$0")"
export NPERM="${NPERM:-200}"
# Use the package's virtual environment (../.venv) unless PYTHON is set explicitly.
if [ -z "${PYTHON:-}" ]; then
  if [ -x ../.venv/bin/python ]; then PYTHON=../.venv/bin/python
  else echo "WARNING: ../.venv not found; falling back to python3 (package versions may differ)" >&2; PYTHON=python3; fi
fi
echo "Using $("$PYTHON" -c 'import sys; print(sys.executable, sys.version.split()[0])')"
for s in 01_series1_descriptors 02_series1_classifiers 03_motif_statistics 06_suppfig2b_ungrouped \
         05_vector_panels 07_motif_panels 09_fig1f_vector 10_suppfig2_composite 08_source_data \
         11_sga_scores 12_sga_panels 17_sga_source_data \
         13_cerevisiae_analysis 14_cerevisiae_montage 15_cerevisiae_source_data \
         16_suppfig4_reproducibility; do
  echo "== $s"; "$PYTHON" "$s.py"
done
echo "All scripts finished."
