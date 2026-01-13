#!/bin/bash
# ========================================================================
# DELFIN Setup Script for BwUniCluster 3.0
# ========================================================================
#
# Automatische Installation von DELFIN auf dem BwUniCluster 3.0
#
# Usage:
#   bash setup_delfin_bwunicluster3.sh
#
# Was macht dieses Skript:
# - Lädt benötigte Module (Python)
# - Klont DELFIN von GitHub (oder aktualisiert bestehendes Repo)
# - Erstellt Virtual Environment
# - Installiert DELFIN im Development-Modus
# - Prüft ORCA-Installation
# - Richtet Shortcuts ein
#
# ========================================================================

set -e  # Bei Fehler abbrechen

# Farben für Output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}"
echo "========================================"
echo "DELFIN Setup für BwUniCluster 3.0"
echo "========================================"
echo -e "${NC}"

# ========================================================================
# Konfiguration - Hier anpassen!
# ========================================================================

PYTHON_MODULE="devel/python/3.11"     # Anpassen an verfügbare Version
ORCA_MODULE="chem/orca/6.1.1"         # Anpassen an verfügbare Version
INSTALL_DIR="$HOME/DELFIN"            # Wo DELFIN installiert wird
VENV_NAME=".venv"                     # Name des Virtual Environments

# ========================================================================
# Funktionen
# ========================================================================

print_step() {
    echo -e "\n${BLUE}▶ $1${NC}"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ $1${NC}"
}

check_command() {
    if command -v $1 &> /dev/null; then
        print_success "$1 gefunden: $(which $1)"
        return 0
    else
        print_error "$1 nicht gefunden"
        return 1
    fi
}

# ========================================================================
# Schritt 1: Module prüfen und laden
# ========================================================================

print_step "Schritt 1: Module prüfen"

echo "Verfügbare Python-Module:"
module avail python 2>&1 | grep -i python || echo "  (module avail python zeigt nichts)"

echo ""
echo "Verfügbare ORCA-Module:"
module avail orca 2>&1 | grep -i orca || echo "  (module avail orca zeigt nichts)"

echo ""
print_step "Lade Python-Modul: $PYTHON_MODULE"
if module load $PYTHON_MODULE 2>/dev/null; then
    print_success "Python-Modul geladen"
    python --version
else
    print_error "Python-Modul '$PYTHON_MODULE' konnte nicht geladen werden"
    echo ""
    echo "Bitte passe PYTHON_MODULE in diesem Skript an:"
    echo "  1. Führe aus: module avail python"
    echo "  2. Finde die richtige Version (z.B. devel/python/3.11)"
    echo "  3. Editiere PYTHON_MODULE=... am Anfang dieses Skripts"
    exit 1
fi

# ========================================================================
# Schritt 2: DELFIN klonen oder aktualisieren
# ========================================================================

print_step "Schritt 2: DELFIN Repository"

if [ -d "$INSTALL_DIR" ]; then
    print_warning "DELFIN bereits vorhanden in: $INSTALL_DIR"
    echo "Möchtest du es aktualisieren? (y/N)"
    read -r -n 1 response
    echo
    if [[ "$response" =~ ^([yY])$ ]]; then
        cd "$INSTALL_DIR"
        print_step "Aktualisiere Repository..."
        git pull
        print_success "Repository aktualisiert"
    else
        print_warning "Überspringe Git-Update"
    fi
else
    print_step "Klone DELFIN Repository..."
    git clone https://github.com/ComPlat/DELFIN.git "$INSTALL_DIR"
    print_success "Repository geklont nach: $INSTALL_DIR"
fi

cd "$INSTALL_DIR"

# ========================================================================
# Schritt 3: Virtual Environment
# ========================================================================

print_step "Schritt 3: Virtual Environment"

if [ -d "$VENV_NAME" ]; then
    print_warning "Virtual Environment existiert bereits"
    echo "Möchtest du es neu erstellen? (y/N)"
    read -r -n 1 response
    echo
    if [[ "$response" =~ ^([yY])$ ]]; then
        rm -rf "$VENV_NAME"
        python -m venv "$VENV_NAME"
        print_success "Virtual Environment neu erstellt"
    else
        print_warning "Verwende bestehendes Virtual Environment"
    fi
else
    python -m venv "$VENV_NAME"
    print_success "Virtual Environment erstellt"
fi

# ========================================================================
# Schritt 4: DELFIN installieren
# ========================================================================

print_step "Schritt 4: DELFIN installieren"

source "$VENV_NAME/bin/activate"
print_success "Virtual Environment aktiviert"

echo "Upgrade pip..."
pip install --upgrade pip -q

echo "Installiere DELFIN im Development-Modus..."
pip install -e . -q

print_success "DELFIN installiert"

# ========================================================================
# Schritt 5: Installation prüfen
# ========================================================================

print_step "Schritt 5: Installation prüfen"

check_command delfin
DELFIN_VERSION=$(delfin --version 2>&1 || echo "unbekannt")
echo -e "${GREEN}DELFIN Version: $DELFIN_VERSION${NC}"

# Git Info
GIT_COMMIT=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")
GIT_BRANCH=$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "unknown")
echo "Git Branch: $GIT_BRANCH"
echo "Git Commit: $GIT_COMMIT"

# ========================================================================
# Schritt 6: ORCA prüfen
# ========================================================================

print_step "Schritt 6: ORCA prüfen"

if module load $ORCA_MODULE 2>/dev/null; then
    print_success "ORCA-Modul geladen: $ORCA_MODULE"
    if check_command orca; then
        ORCA_VERSION=$(orca 2>&1 | grep "Program Version" | head -1 || echo "Version unbekannt")
        echo "$ORCA_VERSION"
    fi
else
    print_warning "ORCA-Modul '$ORCA_MODULE' nicht gefunden"
    echo ""
    echo "ORCA ist erforderlich für DELFIN!"
    echo "Bitte passe ORCA_MODULE in diesem Skript an:"
    echo "  1. Führe aus: module avail orca"
    echo "  2. Finde die richtige Version (z.B. chem/orca/6.1.1)"
    echo "  3. Editiere ORCA_MODULE=... am Anfang dieses Skripts"
fi

# ========================================================================
# Schritt 7: Shortcuts einrichten
# ========================================================================

print_step "Schritt 7: Shortcuts einrichten"

BASHRC="$HOME/.bashrc"
ALIAS_LINE="alias delfin-activate='source $INSTALL_DIR/$VENV_NAME/bin/activate'"
MODULE_LINE="alias delfin-modules='module load $PYTHON_MODULE $ORCA_MODULE'"

if grep -q "delfin-activate" "$BASHRC" 2>/dev/null; then
    print_warning "Alias 'delfin-activate' existiert bereits in ~/.bashrc"
else
    echo "" >> "$BASHRC"
    echo "# DELFIN Shortcuts (automatisch hinzugefügt)" >> "$BASHRC"
    echo "$ALIAS_LINE" >> "$BASHRC"
    echo "$MODULE_LINE" >> "$BASHRC"
    print_success "Shortcuts zu ~/.bashrc hinzugefügt"
fi

# ========================================================================
# Schritt 8: Test-Kommandos
# ========================================================================

print_step "Schritt 8: Erstelle Test-Projekt"

TEST_DIR="$HOME/delfin_test"
if [ ! -d "$TEST_DIR" ]; then
    mkdir -p "$TEST_DIR"
    cd "$TEST_DIR"

    # Erstelle minimales CONTROL.txt
    cat > CONTROL.txt << 'EOF'
# Minimal CONTROL.txt für Test
charge=0
mult=1
method=OCCUPIER
OCCUPIER_method=auto
calc_initial=yes
oxidation_steps=1
reduction_steps=1
pal_jobs=4
EOF

    # Erstelle minimale input.txt (Wassermolekül)
    cat > input.txt << 'EOF'
O     0.0000000000    0.0000000000    0.1173000000
H     0.0000000000    0.7572000000   -0.4692000000
H     0.0000000000   -0.7572000000   -0.4692000000
EOF

    print_success "Test-Projekt erstellt in: $TEST_DIR"
    echo "Zum Testen: cd $TEST_DIR && delfin --help"
fi

# ========================================================================
# Zusammenfassung
# ========================================================================

echo ""
echo -e "${GREEN}"
echo "========================================"
echo "✓ DELFIN Setup erfolgreich!"
echo "========================================"
echo -e "${NC}"

echo -e "\n${BLUE}Installationsdetails:${NC}"
echo "  DELFIN Verzeichnis: $INSTALL_DIR"
echo "  Virtual Environment: $INSTALL_DIR/$VENV_NAME"
echo "  DELFIN Version: $DELFIN_VERSION"
echo "  Git Branch: $GIT_BRANCH"
echo "  Git Commit: $GIT_COMMIT"

echo -e "\n${BLUE}Verwendung:${NC}"
echo "  1. Module laden:"
echo "     ${YELLOW}module load $PYTHON_MODULE $ORCA_MODULE${NC}"
echo ""
echo "  2. DELFIN aktivieren:"
echo "     ${YELLOW}source $INSTALL_DIR/$VENV_NAME/bin/activate${NC}"
echo ""
echo "  3. Oder nutze die Shortcuts (nach neuem Login):"
echo "     ${YELLOW}delfin-modules${NC}      # Lädt Module"
echo "     ${YELLOW}delfin-activate${NC}     # Aktiviert DELFIN"
echo ""
echo "  4. DELFIN verwenden:"
echo "     ${YELLOW}delfin --version${NC}"
echo "     ${YELLOW}delfin --help${NC}"

echo -e "\n${BLUE}Nächste Schritte:${NC}"
echo "  1. Teste die Installation:"
echo "     ${YELLOW}cd $TEST_DIR && delfin-activate && delfin --help${NC}"
echo ""
echo "  2. Erstelle dein Projekt:"
echo "     ${YELLOW}mkdir ~/delfin_projects/mein_projekt${NC}"
echo "     ${YELLOW}cd ~/delfin_projects/mein_projekt${NC}"
echo "     ${YELLOW}delfin --define=struktur.xyz${NC}"
echo ""
echo "  3. Submit-Skript kopieren:"
echo "     ${YELLOW}cp $INSTALL_DIR/examples/example_Job_Submission_Scripts/bwunicluster3_delfin_dev.sh .${NC}"
echo ""
echo "  4. Job starten:"
echo "     ${YELLOW}sbatch bwunicluster3_delfin_dev.sh${NC}"

echo -e "\n${BLUE}Updates durchführen:${NC}"
echo "  ${YELLOW}cd $INSTALL_DIR && git pull${NC}"

echo -e "\n${GREEN}Setup abgeschlossen!${NC}\n"

deactivate 2>/dev/null || true
