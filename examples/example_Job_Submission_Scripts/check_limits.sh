#!/bin/bash
# ========================================================================
# Check BwUniCluster 3.0 Limits
# ========================================================================
#
# Zeigt deine persönlichen Cluster-Limits an
#
# Usage:
#   bash check_limits.sh
#
# ========================================================================

# Farben
BLUE='\033[0;34m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

echo -e "${BLUE}"
echo "========================================"
echo "BwUniCluster 3.0 - Resource Limits"
echo "========================================"
echo -e "${NC}"

# ========================================================================
# 1. Partition Informationen
# ========================================================================

echo -e "\n${CYAN}▶ Verfügbare Partitionen und Limits:${NC}\n"

sinfo -o "%.15P %.10l %.5D %.6c %.8m %.10G %.20C" | head -20

echo -e "\n${YELLOW}Legende:${NC}"
echo "  PARTITION  = Partition-Name"
echo "  TIMELIMIT  = Max. Laufzeit"
echo "  NODES      = Anzahl Nodes in Partition"
echo "  CPUS       = CPUs pro Node"
echo "  MEMORY     = RAM pro Node"
echo "  GRES       = GPUs (falls vorhanden)"
echo "  CPUS(A/I/O/T) = Allocated/Idle/Other/Total"

# ========================================================================
# 2. Deine aktuellen Jobs
# ========================================================================

echo -e "\n${CYAN}▶ Deine laufenden Jobs:${NC}\n"

RUNNING_JOBS=$(squeue -u $USER -t RUNNING -h | wc -l)
PENDING_JOBS=$(squeue -u $USER -t PENDING -h | wc -l)
TOTAL_JOBS=$(squeue -u $USER -h | wc -l)

if [ $TOTAL_JOBS -eq 0 ]; then
    echo "  Keine Jobs aktuell aktiv"
else
    squeue -u $USER -o "%.10i %.15P %.30j %.8u %.2t %.10M %.6D %.5C %.10m %R"

    echo -e "\n${GREEN}Zusammenfassung:${NC}"
    echo "  Laufende Jobs:  $RUNNING_JOBS"
    echo "  Wartende Jobs:  $PENDING_JOBS"
    echo "  Gesamt:         $TOTAL_JOBS"
fi

# Cores berechnen
if [ $RUNNING_JOBS -gt 0 ]; then
    USED_CORES=$(squeue -u $USER -t RUNNING -h -o "%C" | awk '{sum+=$1} END {print sum}')
    echo "  Genutzte Cores: $USED_CORES"
fi

# ========================================================================
# 3. Fair-Share Informationen
# ========================================================================

echo -e "\n${CYAN}▶ Deine Fair-Share Priorität:${NC}\n"

if command -v sshare &> /dev/null; then
    sshare -u $USER -o "Account,User,RawShares,NormShares,RawUsage,EffectvUsage,FairShare" 2>/dev/null || echo "  Fair-Share Info nicht verfügbar"

    echo -e "\n${YELLOW}Fair-Share Bedeutung:${NC}"
    echo "  FairShare > 0.5  = Hohe Priorität (du kannst mehr nutzen)"
    echo "  FairShare ≈ 0.5  = Normale Priorität"
    echo "  FairShare < 0.5  = Niedrige Priorität (du hast viel genutzt)"
else
    echo "  sshare nicht verfügbar"
fi

# ========================================================================
# 4. Account Informationen
# ========================================================================

echo -e "\n${CYAN}▶ Dein Account:${NC}\n"

# Sacctmgr für detaillierte Limits
if command -v sacctmgr &> /dev/null; then
    echo "Account Details:"
    sacctmgr show user $USER format=User,Account,Partition,MaxJobs,MaxSubmit,MaxCPUs,MaxNodes -p 2>/dev/null | column -t -s '|' || echo "  Account-Limits nicht verfügbar"
else
    echo "  sacctmgr nicht verfügbar"
fi

# ========================================================================
# 5. Partition-spezifische Limits
# ========================================================================

echo -e "\n${CYAN}▶ Detaillierte Partition-Limits:${NC}\n"

for partition in cpu cpu_il highmem dev_cpu; do
    echo -e "${YELLOW}Partition: $partition${NC}"
    scontrol show partition $partition 2>/dev/null | grep -E "(MaxTime|MaxNodes|MaxCPUsPerUser|DefMemPerCPU|MaxMemPerNode)" || echo "  Partition nicht gefunden"
    echo ""
done

# ========================================================================
# 6. QOS (Quality of Service) Limits
# ========================================================================

echo -e "\n${CYAN}▶ Quality of Service (QOS) Limits:${NC}\n"

if command -v sacctmgr &> /dev/null; then
    sacctmgr show qos format=Name,MaxWall,MaxJobsPerUser,MaxSubmitJobsPerUser,MaxCPUsPerUser -p 2>/dev/null | column -t -s '|' || echo "  QOS-Limits nicht verfügbar"
else
    echo "  QOS-Info nicht verfügbar"
fi

# ========================================================================
# 7. Praktische Empfehlungen
# ========================================================================

echo -e "\n${CYAN}▶ Praktische Empfehlungen für DELFIN:${NC}\n"

# Cluster-Auslastung
TOTAL_IDLE_CORES=$(sinfo -h -p cpu -o "%C" | cut -d'/' -f2)
TOTAL_CORES=$(sinfo -h -p cpu -o "%C" | cut -d'/' -f4)

if [ ! -z "$TOTAL_IDLE_CORES" ] && [ ! -z "$TOTAL_CORES" ]; then
    USAGE_PERCENT=$((100 * (TOTAL_CORES - TOTAL_IDLE_CORES) / TOTAL_CORES))
    echo "Cluster-Auslastung (cpu partition): ${USAGE_PERCENT}%"

    if [ $USAGE_PERCENT -lt 50 ]; then
        echo -e "${GREEN}✓ Cluster wenig ausgelastet - gute Zeit für große Jobs${NC}"
        echo "  Empfehlung: 10-20 Jobs à 40 Cores (400-800 Cores)"
    elif [ $USAGE_PERCENT -lt 80 ]; then
        echo -e "${YELLOW}⚠ Cluster moderat ausgelastet${NC}"
        echo "  Empfehlung: 5-10 Jobs à 40 Cores (200-400 Cores)"
    else
        echo -e "${YELLOW}⚠ Cluster stark ausgelastet${NC}"
        echo "  Empfehlung: 3-5 Jobs à 40 Cores (120-200 Cores)"
    fi
else
    echo "Auslastung konnte nicht ermittelt werden"
fi

echo -e "\n${YELLOW}Allgemeine Empfehlungen:${NC}"
echo "  - Start: 3-5 Jobs (120-200 Cores)"
echo "  - Standard: 10-15 Jobs (400-600 Cores)"
echo "  - Maximum: 20-30 Jobs (800-1200 Cores)"
echo "  - dev_cpu: Nur für Tests, max 1 Job, 30 Min"

# ========================================================================
# 8. Wie viele Jobs kannst du noch starten?
# ========================================================================

echo -e "\n${CYAN}▶ Freie Ressourcen (cpu partition):${NC}\n"

# Freie Nodes
FREE_NODES=$(sinfo -h -p cpu -o "%a/%D" | head -1)
echo "Nodes: $FREE_NODES (verfügbar/gesamt)"

# Freie CPUs
CPU_INFO=$(sinfo -h -p cpu -o "%C" | head -1)
echo "CPUs: $CPU_INFO (allocated/idle/other/total)"

# Wartezeit schätzen
if [ $PENDING_JOBS -gt 0 ]; then
    echo -e "\n${YELLOW}Du hast $PENDING_JOBS wartende Jobs${NC}"
    echo "Geschätzte Wartezeit: $(squeue -u $USER -t PENDING -o "%.10i %.20S" 2>/dev/null | tail -n +2 | head -1 | awk '{print $2}')"
fi

# ========================================================================
# 9. Nützliche Kommandos
# ========================================================================

echo -e "\n${CYAN}▶ Nützliche Kommandos:${NC}\n"

cat << 'EOF'
# Deine Jobs anzeigen
squeue -u $USER

# Partition-Auslastung
sinfo -p cpu

# Detaillierte Job-Info
scontrol show job <job_id>

# Fair-Share Priorität
sshare -u $USER

# Account-Limits
sacctmgr show user $USER

# Geschätzte Startzeit
squeue -u $USER --start

# Cluster-Statistiken
squeue --start | head -20
EOF

echo ""
echo -e "${GREEN}========================================"
echo "Limit-Check abgeschlossen!"
echo "========================================${NC}"
echo ""
