#!/usr/bin/env bash
# =============================================================================
# run_all.sh -- run the R pipeline end to end on the simulated ground-truth
# dataset (or on real data with --raw once data/raw/ is populated).
#
#   ./run_all.sh                 # simulate, then run every step
#   ./run_all.sh --raw           # skip simulation, use data/raw/
#   ./run_all.sh --config config/config_prior.yaml
#
# On Windows the Rscript path is auto-detected; override with RSCRIPT=...
# =============================================================================
set -euo pipefail

CONFIG="config/config.yaml"
USE_RAW=""
while [ $# -gt 0 ]; do
  case "$1" in
    --raw)    USE_RAW="--raw"; shift ;;
    --config) CONFIG="$2"; shift 2 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

if [ -z "${RSCRIPT:-}" ]; then
  if command -v Rscript >/dev/null 2>&1; then
    RSCRIPT="Rscript"
  else
    RSCRIPT="$(ls -d /c/Program\ Files/R/R-*/bin/Rscript.exe 2>/dev/null | tail -1 || true)"
  fi
fi
[ -n "$RSCRIPT" ] || { echo "Rscript not found; set RSCRIPT=/path/to/Rscript" >&2; exit 1; }
echo "Rscript : $RSCRIPT"
echo "config  : $CONFIG"

step () { echo; echo "=============== $1 ==============="; shift; "$RSCRIPT" "$@"; }

if [ -z "$USE_RAW" ]; then
  step "00  simulate Case B data"      R/00_simulate_data.R      "$CONFIG"
  step "00b validate the simulator"    R/00b_check_simulation.R  "$CONFIG"
fi
step "01  QC and normalisation"        R/01_qc_normalise.R       "$CONFIG" $USE_RAW
step "02  ID mapping"                  R/02_id_mapping.R         "$CONFIG"
step "03  missingness"                 R/03_missingness.R        "$CONFIG"
step "04  univariate temporal DE"      R/04_univariate_temporal_de.R "$CONFIG"
step "05  concordance and archetypes"  R/05_concordance_archetypes.R "$CONFIG"
step "06  Case B integration"          R/06_integration_caseB.R  "$CONFIG"

echo
echo "=============== done ==============="
echo "tables : results/{qc,mapping,de,concordance,integration}/"
echo "figures: figures/"
