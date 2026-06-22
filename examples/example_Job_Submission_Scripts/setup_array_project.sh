#!/bin/bash
# ========================================================================
# Setup DELFIN Job Array Project
# ========================================================================
#
# Erstellt automatisch die Verzeichnisstruktur für Job Arrays
#
# Usage:
#   bash setup_array_project.sh <anzahl_moleküle>
#
# Beispiel:
#   bash setup_array_project.sh 10
#
# ========================================================================

# Farben
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Anzahl Moleküle
NUM_MOLECULES=${1:-10}  # Default: 10

echo -e "${BLUE}"
echo "========================================"
echo "DELFIN Array Project Setup"
echo "========================================"
echo -e "${NC}"

echo "Erstelle Struktur für $NUM_MOLECULES Moleküle..."
echo ""

# Hauptverzeichnis erstellen
PROJECT_DIR="delfin_array_project"
if [ -d "$PROJECT_DIR" ]; then
    echo -e "${YELLOW}Warnung: $PROJECT_DIR existiert bereits${NC}"
    read -p "Fortfahren? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Abgebrochen."
        exit 1
    fi
else
    mkdir -p "$PROJECT_DIR"
fi

cd "$PROJECT_DIR"

# Molekül-Verzeichnisse erstellen
for i in $(seq 1 $NUM_MOLECULES); do
    MOL_DIR="molecule_$i"

    if [ ! -d "$MOL_DIR" ]; then
        mkdir -p "$MOL_DIR"

        # Minimal CONTROL.txt
        cat > "$MOL_DIR/CONTROL.txt" << 'EOF'
# CONTROL.txt für Molekül (anpassen!)
charge=0
mult=1
input_file=input.txt
method=OCCUPIER
OCCUPIER_method=auto
OCCUPIER_tree=deep
calc_initial=yes
oxidation_steps=1,2
reduction_steps=1,2
parallel_workflows=yes
enable_auto_recovery=yes
max_recovery_attempts=1
enable_job_timeouts=no
XTB_OPT=yes
method_name=B3LYP
dispersion=D3BJ
basis_set=def2-SVP
basis_set_FSPE=def2-TZVP
EOF

        # Minimal input.txt (Wasser als Platzhalter)
        cat > "$MOL_DIR/input.txt" << 'EOF'
O     0.0000000000    0.0000000000    0.1173000000
H     0.0000000000    0.7572000000   -0.4692000000
H     0.0000000000   -0.7572000000   -0.4692000000
EOF

        echo -e "${GREEN}✓${NC} $MOL_DIR erstellt"
    else
        echo -e "${YELLOW}⚠${NC} $MOL_DIR existiert bereits (übersprungen)"
    fi
done

echo ""
echo -e "${GREEN}Projekt-Struktur erstellt!${NC}"
echo ""
echo -e "${BLUE}Verzeichnis-Struktur:${NC}"
echo "$PROJECT_DIR/"
echo "  ├── molecule_1/"
echo "  │   ├── CONTROL.txt"
echo "  │   └── input.txt"
echo "  ├── molecule_2/"
echo "  │   ├── CONTROL.txt"
echo "  │   └── input.txt"
echo "  ..."
echo "  └── molecule_$NUM_MOLECULES/"
echo "      ├── CONTROL.txt"
echo "      └── input.txt"
echo ""

echo -e "${BLUE}Nächste Schritte:${NC}"
echo ""
echo "1. Ersetze die Platzhalter-Geometrien:"
echo "   ${YELLOW}cd $PROJECT_DIR/molecule_1${NC}"
echo "   ${YELLOW}# Füge deine XYZ-Koordinaten in input.txt ein${NC}"
echo "   ${YELLOW}nano input.txt${NC}"
echo ""
echo "2. Passe CONTROL.txt an (charge, mult, etc.):"
echo "   ${YELLOW}nano CONTROL.txt${NC}"
echo ""
echo "3. Wiederhole für alle Moleküle (molecule_2, molecule_3, ...)"
echo ""
echo "4. Kopiere Submit-Skript:"
echo "   ${YELLOW}cp ~/UniBW/submit_array_small.sh $PROJECT_DIR/${NC}"
echo ""
echo "5. Passe Array-Größe im Skript an:"
echo "   ${YELLOW}nano submit_array_small.sh${NC}"
echo "   ${YELLOW}# Ändere: #SBATCH --array=1-$NUM_MOLECULES%5${NC}"
echo ""
echo "6. Starte Job Array:"
echo "   ${YELLOW}cd $PROJECT_DIR${NC}"
echo "   ${YELLOW}sbatch submit_array_small.sh${NC}"
echo ""
echo "7. Überwache Jobs:"
echo "   ${YELLOW}squeue -u \$USER${NC}"
echo "   ${YELLOW}tail -f delfin_*.out${NC}"
echo ""

echo -e "${GREEN}Setup abgeschlossen!${NC}"
