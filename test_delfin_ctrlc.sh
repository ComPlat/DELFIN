#!/bin/bash
# Test Ctrl+C behavior with DELFIN

cd /home/qmchem_max/ComPlat/TEST_DELFIN/DELFIN_OCCUPIER_p_own2

# Clean up
rm -rf initial_OCCUPIER/ *.log *.out *.inp 2>/dev/null

echo "Starting DELFIN..."
echo "Wait for 'INFO: Starting ORCA...' then press Ctrl+C"
echo "You should see '🛑 SIGINT received' within 1 second"
echo ""

delfin
