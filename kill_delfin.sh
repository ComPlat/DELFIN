#!/bin/bash
# Helper script to manage DELFIN processes

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}=== DELFIN Process Manager ===${NC}\n"

# Function to get working directory of a process
get_cwd() {
    local pid=$1
    if [ -d "/proc/$pid" ]; then
        readlink -f "/proc/$pid/cwd" 2>/dev/null || echo "unknown"
    else
        echo "not running"
    fi
}

# Find all DELFIN processes
mapfile -t PIDS < <(pgrep -f "python.*delfin" 2>/dev/null || true)

if [ ${#PIDS[@]} -eq 0 ]; then
    echo -e "${GREEN}No DELFIN processes running${NC}"
    exit 0
fi

echo -e "${YELLOW}Found ${#PIDS[@]} DELFIN process(es):${NC}\n"

# Show all processes with their working directories
for i in "${!PIDS[@]}"; do
    pid="${PIDS[$i]}"
    cwd=$(get_cwd "$pid")

    # Get process start time
    start_time=$(ps -p "$pid" -o lstart= 2>/dev/null || echo "unknown")

    # Get ORCA child processes
    orca_count=$(pgrep -P "$pid" | xargs -I {} pgrep -P {} 2>/dev/null | wc -l || echo "0")

    echo -e "${GREEN}[$((i+1))]${NC} PID: ${BLUE}$pid${NC}"
    echo "    Directory: $cwd"
    echo "    Started: $start_time"
    echo "    ORCA children: $orca_count"
    echo ""
done

# Ask user which to kill
if [ "$1" == "--kill-all" ]; then
    echo -e "${RED}Killing ALL DELFIN processes...${NC}"
    for pid in "${PIDS[@]}"; do
        echo "Killing PID $pid..."
        kill -9 "$pid" 2>/dev/null || true
    done
    echo -e "${GREEN}Done!${NC}"
    exit 0
fi

if [ "$1" == "--kill" ]; then
    if [ -z "$2" ]; then
        echo -e "${RED}Error: Please specify process number (1-${#PIDS[@]})${NC}"
        echo "Usage: $0 --kill <number>"
        exit 1
    fi

    idx=$((${2} - 1))
    if [ $idx -lt 0 ] || [ $idx -ge ${#PIDS[@]} ]; then
        echo -e "${RED}Error: Invalid process number${NC}"
        exit 1
    fi

    pid="${PIDS[$idx]}"
    echo -e "${YELLOW}Killing DELFIN process $pid and all its ORCA children...${NC}"

    # Kill all child processes (including ORCA)
    pkill -9 -P "$pid" 2>/dev/null || true

    # Kill the DELFIN process itself
    kill -9 "$pid" 2>/dev/null || true

    echo -e "${GREEN}Process $pid killed!${NC}"
    exit 0
fi

if [ "$1" == "--kill-dir" ]; then
    if [ -z "$2" ]; then
        echo -e "${RED}Error: Please specify directory path${NC}"
        echo "Usage: $0 --kill-dir <directory>"
        exit 1
    fi

    target_dir=$(readlink -f "$2")
    found=0

    for pid in "${PIDS[@]}"; do
        cwd=$(get_cwd "$pid")
        if [ "$cwd" == "$target_dir" ]; then
            echo -e "${YELLOW}Killing DELFIN in $target_dir (PID $pid)...${NC}"
            pkill -9 -P "$pid" 2>/dev/null || true
            kill -9 "$pid" 2>/dev/null || true
            echo -e "${GREEN}Killed!${NC}"
            found=1
        fi
    done

    if [ $found -eq 0 ]; then
        echo -e "${RED}No DELFIN process found in $target_dir${NC}"
        exit 1
    fi
    exit 0
fi

echo -e "${YELLOW}Usage:${NC}"
echo "  $0                      # Show all DELFIN processes"
echo "  $0 --kill <number>      # Kill specific process (1-${#PIDS[@]})"
echo "  $0 --kill-dir <path>    # Kill DELFIN running in specific directory"
echo "  $0 --kill-all           # Kill ALL DELFIN processes"
