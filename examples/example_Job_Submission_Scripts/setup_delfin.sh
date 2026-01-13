#!/bin/bash
# ========================================================================
# DELFIN Setup Script for BwUniCluster 3.0
# ========================================================================
#
# Automatische Installation von DELFIN auf dem BwUniCluster 3.0
# - Erkennt verfügbare Module automatisch
# - Interaktive Auswahl von Python und ORCA
# - Passt alle Submit-Skripte automatisch an
#
# Usage:
#   bash setup_delfin.sh
#
# ========================================================================

set -e  # Bei Fehler abbrechen

# Farben für Output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

echo -e "${BLUE}"
echo "========================================"
echo "DELFIN Setup für BwUniCluster 3.0"
echo "========================================"
echo -e "${NC}"

# ========================================================================
# Konfiguration
# ========================================================================

INSTALL_DIR="$HOME/DELFIN"
VENV_NAME=".venv"
UNIBW_DIR="$HOME/UniBW"
CONFIG_FILE="$UNIBW_DIR/.delfin_config"

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

print_info() {
    echo -e "${CYAN}ℹ $1${NC}"
}

# ========================================================================
# Schritt 1: Module automatisch erkennen
# ========================================================================

print_step "Schritt 1: Verfügbare Module erkennen"

print_info "Suche Python-Module..."
PYTHON_MODULES=($(module avail python 2>&1 | grep -Eo '(python|devel/python|lang/python)[^ ]*' | sort -V))

if [ ${#PYTHON_MODULES[@]} -eq 0 ]; then
    print_error "Keine Python-Module gefunden!"
    echo "Bitte prüfe manuell: module avail python"
    exit 1
fi

print_info "Suche ORCA-Module..."
ORCA_MODULES=($(module avail orca 2>&1 | grep -Eo '(orca|chem/orca)[^ ]*' | sort -V))

if [ ${#ORCA_MODULES[@]} -eq 0 ]; then
    print_error "Keine ORCA-Module gefunden!"
    echo "Bitte prüfe manuell: module avail orca"
    exit 1
fi

# ========================================================================
# Schritt 2: Interaktive Modul-Auswahl
# ========================================================================

print_step "Schritt 2: Modul-Auswahl"

echo -e "\n${CYAN}Verfügbare Python-Module:${NC}"
for i in "${!PYTHON_MODULES[@]}"; do
    echo "  [$i] ${PYTHON_MODULES[$i]}"
done

echo -e "\n${YELLOW}Welches Python-Modul möchtest du verwenden?${NC}"
if [ ${#PYTHON_MODULES[@]} -eq 1 ]; then
    PYTHON_CHOICE=0
    echo "Automatisch gewählt: ${PYTHON_MODULES[0]}"
else
    # Vorschlag: letzte Version
    DEFAULT_PYTHON=$((${#PYTHON_MODULES[@]} - 1))
    read -p "Nummer eingeben [Standard: $DEFAULT_PYTHON]: " PYTHON_CHOICE
    PYTHON_CHOICE=${PYTHON_CHOICE:-$DEFAULT_PYTHON}
fi

PYTHON_MODULE="${PYTHON_MODULES[$PYTHON_CHOICE]}"
print_success "Python-Modul gewählt: $PYTHON_MODULE"

echo -e "\n${CYAN}Verfügbare ORCA-Module:${NC}"
for i in "${!ORCA_MODULES[@]}"; do
    echo "  [$i] ${ORCA_MODULES[$i]}"
done

echo -e "\n${YELLOW}Welches ORCA-Modul möchtest du verwenden?${NC}"
if [ ${#ORCA_MODULES[@]} -eq 1 ]; then
    ORCA_CHOICE=0
    echo "Automatisch gewählt: ${ORCA_MODULES[0]}"
else
    # Vorschlag: letzte Version
    DEFAULT_ORCA=$((${#ORCA_MODULES[@]} - 1))
    read -p "Nummer eingeben [Standard: $DEFAULT_ORCA]: " ORCA_CHOICE
    ORCA_CHOICE=${ORCA_CHOICE:-$DEFAULT_ORCA}
fi

ORCA_MODULE="${ORCA_MODULES[$ORCA_CHOICE]}"
print_success "ORCA-Modul gewählt: $ORCA_MODULE"

# Konfiguration speichern
echo "PYTHON_MODULE=$PYTHON_MODULE" > "$CONFIG_FILE"
echo "ORCA_MODULE=$ORCA_MODULE" >> "$CONFIG_FILE"
print_success "Konfiguration gespeichert: $CONFIG_FILE"

# ========================================================================
# Schritt 3: Module laden und testen
# ========================================================================

print_step "Schritt 3: Module testen"

if module load $PYTHON_MODULE 2>/dev/null; then
    print_success "Python-Modul geladen: $PYTHON_MODULE"
    PYTHON_VERSION=$(python --version 2>&1)
    print_info "$PYTHON_VERSION"
else
    print_error "Python-Modul konnte nicht geladen werden"
    exit 1
fi

if module load $ORCA_MODULE 2>/dev/null; then
    print_success "ORCA-Modul geladen: $ORCA_MODULE"
    if command -v orca &> /dev/null; then
        ORCA_PATH=$(which orca)
        print_info "ORCA gefunden: $ORCA_PATH"
    else
        print_warning "ORCA nicht im PATH, aber Modul geladen"
    fi
else
    print_warning "ORCA-Modul konnte nicht geladen werden"
    print_warning "ORCA muss für DELFIN verfügbar sein!"
fi

# ========================================================================
# Schritt 4: DELFIN klonen oder aktualisieren
# ========================================================================

print_step "Schritt 4: DELFIN Repository"

if [ -d "$INSTALL_DIR" ]; then
    print_warning "DELFIN bereits vorhanden: $INSTALL_DIR"
    read -p "Aktualisieren? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        cd "$INSTALL_DIR"
        git pull
        print_success "Repository aktualisiert"
    fi
else
    git clone https://github.com/ComPlat/DELFIN.git "$INSTALL_DIR"
    print_success "Repository geklont: $INSTALL_DIR"
fi

cd "$INSTALL_DIR"

# ========================================================================
# Schritt 5: Virtual Environment
# ========================================================================

print_step "Schritt 5: Virtual Environment"

if [ -d "$VENV_NAME" ]; then
    print_warning "Virtual Environment existiert bereits"
    read -p "Neu erstellen? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        rm -rf "$VENV_NAME"
        python -m venv "$VENV_NAME"
        print_success "Virtual Environment neu erstellt"
    fi
else
    python -m venv "$VENV_NAME"
    print_success "Virtual Environment erstellt"
fi

# ========================================================================
# Schritt 6: DELFIN installieren
# ========================================================================

print_step "Schritt 6: DELFIN installieren"

source "$VENV_NAME/bin/activate"
print_info "Upgrade pip..."
pip install --upgrade pip -q
print_info "Installiere DELFIN..."
pip install -e . -q
print_success "DELFIN installiert"

# ========================================================================
# Schritt 7: Installation prüfen
# ========================================================================

print_step "Schritt 7: Installation prüfen"

DELFIN_VERSION=$(delfin --version 2>&1 || echo "unbekannt")
echo -e "${GREEN}DELFIN Version: $DELFIN_VERSION${NC}"
echo "Git Branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Git Commit: $(git rev-parse --short HEAD)"

# ========================================================================
# Schritt 8: Submit-Skripte automatisch anpassen
# ========================================================================

print_step "Schritt 8: Submit-Skripte anpassen"

if [ -d "$UNIBW_DIR" ]; then
    SUBMIT_SCRIPTS=(
        "$UNIBW_DIR/submit_standard.sh"
        "$UNIBW_DIR/submit_highmem.sh"
        "$UNIBW_DIR/submit_test.sh"
    )

    for script in "${SUBMIT_SCRIPTS[@]}"; do
        if [ -f "$script" ]; then
            # Backup erstellen
            cp "$script" "${script}.bak"

            # Module-Namen ersetzen
            sed -i "s|module load chem/orca/[^ ]*|module load $ORCA_MODULE|g" "$script"
            sed -i "s|module load devel/python/[^ ]*|module load $PYTHON_MODULE|g" "$script"
            sed -i "s|module load python/[^ ]*|module load $PYTHON_MODULE|g" "$script"

            print_success "$(basename $script) angepasst"
        fi
    done

    print_info "Backups wurden erstellt (*.bak)"
else
    print_warning "UniBW-Verzeichnis nicht gefunden: $UNIBW_DIR"
fi

# ========================================================================
# Schritt 9: Shortcuts einrichten
# ========================================================================

print_step "Schritt 9: Shortcuts einrichten"

BASHRC="$HOME/.bashrc"
ALIAS1="alias delfin-activate='source $INSTALL_DIR/$VENV_NAME/bin/activate'"
ALIAS2="alias delfin-modules='module load $PYTHON_MODULE $ORCA_MODULE'"

if grep -q "delfin-activate" "$BASHRC" 2>/dev/null; then
    print_warning "Shortcuts existieren bereits in ~/.bashrc"
    read -p "Aktualisieren? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        # Alte Einträge entfernen
        sed -i '/# DELFIN Shortcuts/d' "$BASHRC"
        sed -i '/delfin-activate/d' "$BASHRC"
        sed -i '/delfin-modules/d' "$BASHRC"

        # Neue Einträge hinzufügen
        echo "" >> "$BASHRC"
        echo "# DELFIN Shortcuts" >> "$BASHRC"
        echo "$ALIAS1" >> "$BASHRC"
        echo "$ALIAS2" >> "$BASHRC"
        print_success "Shortcuts aktualisiert"
    fi
else
    echo "" >> "$BASHRC"
    echo "# DELFIN Shortcuts" >> "$BASHRC"
    echo "$ALIAS1" >> "$BASHRC"
    echo "$ALIAS2" >> "$BASHRC"
    print_success "Shortcuts zu ~/.bashrc hinzugefügt"
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
echo "  DELFIN:         $INSTALL_DIR"
echo "  Python-Modul:   $PYTHON_MODULE"
echo "  ORCA-Modul:     $ORCA_MODULE"
echo "  DELFIN Version: $DELFIN_VERSION"
echo "  Konfiguration:  $CONFIG_FILE"

echo -e "\n${BLUE}Verwendung:${NC}"
echo "  ${YELLOW}delfin-modules${NC}      # Lädt: $PYTHON_MODULE $ORCA_MODULE"
echo "  ${YELLOW}delfin-activate${NC}     # Aktiviert DELFIN Environment"
echo "  ${YELLOW}delfin --version${NC}    # Zeigt DELFIN Version"

echo -e "\n${BLUE}Angepasste Submit-Skripte:${NC}"
echo "  ${GREEN}✓${NC} $UNIBW_DIR/submit_standard.sh"
echo "  ${GREEN}✓${NC} $UNIBW_DIR/submit_highmem.sh"
echo "  ${GREEN}✓${NC} $UNIBW_DIR/submit_test.sh"
echo -e "  ${CYAN}(Backups: *.bak)${NC}"

echo -e "\n${BLUE}Nächste Schritte:${NC}"
echo "  1. Neue Shell öffnen oder: ${YELLOW}source ~/.bashrc${NC}"
echo ""
echo "  2. Projekt erstellen:"
echo "     ${YELLOW}mkdir ~/delfin_projects/mein_molekuel${NC}"
echo "     ${YELLOW}cd ~/delfin_projects/mein_molekuel${NC}"
echo ""
echo "  3. DELFIN einrichten:"
echo "     ${YELLOW}delfin-modules${NC}"
echo "     ${YELLOW}delfin-activate${NC}"
echo "     ${YELLOW}delfin --define=struktur.xyz${NC}"
echo ""
echo "  4. Submit-Skript kopieren:"
echo "     ${YELLOW}cp $UNIBW_DIR/submit_standard.sh .${NC}"
echo ""
echo "  5. Job starten:"
echo "     ${YELLOW}sbatch submit_standard.sh${NC}"

echo -e "\n${BLUE}Updates durchführen:${NC}"
echo "  ${YELLOW}cd $INSTALL_DIR && git pull${NC}"

echo -e "\n${GREEN}Setup abgeschlossen! 🚀${NC}\n"

deactivate 2>/dev/null || true
