#!/bin/bash
# ========================================================================
# DELFIN Local Run Script (no Slurm)
# ========================================================================
#
# MODES (set via DELFIN_MODE environment variable):
#   delfin | delfin-recalc | delfin-recalc-override | orca | build | delfin-co2-chain | auto
#
# ENVIRONMENT VARIABLES:
#   DELFIN_MODE          : Run mode (default: auto)
#   DELFIN_JOB_NAME      : Display name
#   DELFIN_INP_FILE      : Specific .inp file for ORCA mode
#   DELFIN_OVERRIDE      : Override value for delfin-recalc-override
#   DELFIN_JOB_ID        : Job ID for status tracking
#   BUILD_MULTIPLICITY   : Spin multiplicity for build mode (default: 1)
#   DELFIN_CO2_SPECIES_DELTA : Redox species delta for delfin-co2-chain mode (default: 0)
#
# ========================================================================

set -euo pipefail

# === Configuration from environment variables ===
MODE="${DELFIN_MODE:-auto}"
INP_FILE="${DELFIN_INP_FILE:-}"
OVERRIDE="${DELFIN_OVERRIDE:-}"
DISPLAY_JOB_NAME="${DELFIN_JOB_NAME:-local_job}"
JOB_ID="${DELFIN_JOB_ID:-0}"

# === MPI Configuration for ORCA stability ===
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,tcp,vader
export OMPI_MCA_mpi_show_mca_params_file=0
export OMPI_MCA_mpi_yield_when_idle=1
export OMPI_MCA_coll_hcoll_enable=0
unset OMPI_MCA_btl_tcp_if_include 2>/dev/null || true
unset OMPI_MCA_btl_tcp_if_exclude 2>/dev/null || true
export OMPI_MCA_hwloc_base_binding_policy=none
export OMPI_MCA_rmaps_base_mapping_policy=core
export OMPI_MCA_rmaps_base_oversubscribe=true

# === Environment ===
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export DELFIN_ORCA_PROGRESS=0
export MPLBACKEND=Agg

# Job info
echo "========================================"
echo "DELFIN Local Job - $DISPLAY_JOB_NAME"
echo "Mode: $MODE"
echo "========================================"
echo "Job ID:      $JOB_ID"
echo "Host:        $(hostname)"
echo "CPUs avail:  $(nproc)"
echo "Working Dir: $PWD"
echo "Started:     $(date)"
if [ -n "$INP_FILE" ]; then
    echo "Input File:  $INP_FILE"
fi
echo "========================================"
echo ""

# Check tools
echo "ORCA: $(command -v orca || echo 'not found')"
echo "mpirun: $(command -v mpirun || echo 'not found')"
echo "delfin: $(command -v delfin || echo 'not found')"
echo ""

# Auto-detect mode if set to "auto"
if [ "$MODE" = "auto" ]; then
    if [ -f "CONTROL.txt" ] && [ -f "input.txt" ]; then
        MODE="delfin"
        echo "Auto-detected mode: DELFIN (CONTROL.txt + input.txt found)"
    elif ls *.inp 1>/dev/null 2>&1; then
        MODE="orca"
        echo "Auto-detected mode: ORCA (*.inp files found)"
    else
        echo "ERROR: No valid input files found"
        echo "       Expected: CONTROL.txt + input.txt (DELFIN) or *.inp (ORCA)"
        echo 1 > ".exit_code_${JOB_ID}"
        exit 1
    fi
fi

# Run appropriate mode
EXIT_CODE=0
set +e
case "$MODE" in
    delfin)
        echo "Starting DELFIN..."
        delfin
        EXIT_CODE=$?
        ;;
    delfin-recalc)
        echo "Starting DELFIN --recalc..."
        delfin --recalc
        EXIT_CODE=$?
        ;;
    delfin-recalc-override)
        if [ -z "$OVERRIDE" ]; then
            echo "ERROR: DELFIN_OVERRIDE not set for delfin-recalc-override mode"
            EXIT_CODE=1
        else
            echo "Starting DELFIN --recalc --occupier-override $OVERRIDE..."
            delfin --recalc --occupier-override "$OVERRIDE"
            EXIT_CODE=$?
        fi
        ;;
    orca)
        if [ -z "$INP_FILE" ]; then
            INP_FILE=$(ls *.inp 2>/dev/null | head -1)
        fi
        if [ -z "$INP_FILE" ]; then
            echo "ERROR: No .inp file found for ORCA mode"
            EXIT_CODE=1
        elif [ ! -f "$INP_FILE" ]; then
            echo "ERROR: Input file not found: $INP_FILE"
            EXIT_CODE=1
        else
            OUT_FILE="${INP_FILE%.inp}.out"
            echo "Starting ORCA: $INP_FILE -> $OUT_FILE"
            orca "$INP_FILE" > "$OUT_FILE" 2>&1
            EXIT_CODE=$?
        fi
        ;;
    build)
        BUILD_MULT="${BUILD_MULTIPLICITY:-1}"
        echo "Starting delfin-build (complex build-up)..."
        echo "  Multiplicity: $BUILD_MULT"
        "$DELFIN_PYTHON" -m delfin.build_up_complex input.txt --goat --directory builder --multiplicity "$BUILD_MULT" --verbose
        EXIT_CODE=$?
        ;;
    delfin-co2-chain)
        CO2_SPECIES_DELTA="${DELFIN_CO2_SPECIES_DELTA:-0}"
        echo "Starting DELFIN + CO2 Coordinator chain..."
        echo "  CO2 Species Delta: $CO2_SPECIES_DELTA"

        # Resolve the Python that belongs to the delfin installation
        # (bare `python` may point to a different env like chemdarwin)
        DELFIN_PYTHON="$(head -1 "$(command -v delfin)" | sed 's/^#!//')"
        if [ ! -x "$DELFIN_PYTHON" ]; then
            DELFIN_PYTHON="$(command -v python)"
        fi

        # Step 1: DELFIN
        echo "=== Step 1/2: Running DELFIN ==="
        delfin
        DELFIN_EXIT=$?
        if [ "$DELFIN_EXIT" -ne 0 ]; then
            echo "ERROR: DELFIN failed (exit $DELFIN_EXIT). CO2 Coordinator skipped."
            EXIT_CODE=$DELFIN_EXIT
        else
            # Step 2: Chain Setup + CO2 Coordinator
            echo "=== Step 2/2: CO2 Coordinator ==="
            "$DELFIN_PYTHON" -m delfin.co2.chain_setup "$CO2_SPECIES_DELTA"
            SETUP_EXIT=$?
            if [ "$SETUP_EXIT" -ne 0 ]; then
                echo "ERROR: CO2 setup failed (exit $SETUP_EXIT)"
                EXIT_CODE=$SETUP_EXIT
            else
                cd CO2_coordination
                delfin co2
                EXIT_CODE=$?
                cd ..
            fi
        fi
        ;;
    *)
        echo "ERROR: Unknown mode: $MODE"
        EXIT_CODE=1
        ;;
esac
set -e

echo ""
echo "========================================"
echo "Job finished: $(date)"
echo "Exit Code:   $EXIT_CODE"
echo "========================================"

# Write exit code file for status tracking
echo "$EXIT_CODE" > ".exit_code_${JOB_ID}"

exit $EXIT_CODE
