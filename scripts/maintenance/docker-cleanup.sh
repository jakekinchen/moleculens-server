#!/usr/bin/env bash
# Docker Cleanup Script - Prevents disk space issues
# Usage: ./docker-cleanup.sh [--force] [--threshold=80]

set -e

# Default settings
FORCE_CLEANUP=false
DISK_THRESHOLD=80
DRY_RUN=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --force)
            FORCE_CLEANUP=true
            shift
            ;;
        --threshold=*)
            DISK_THRESHOLD="${1#*=}"
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [--force] [--threshold=80] [--dry-run]"
            echo "  --force      : Clean up regardless of disk usage"
            echo "  --threshold  : Disk usage percentage to trigger cleanup (default: 80)"
            echo "  --dry-run    : Show what would be cleaned without doing it"
            exit 0
            ;;
        *)
            echo "Unknown option $1"
            exit 1
            ;;
    esac
done

# Function to get disk usage percentage
get_disk_usage() {
    df / | awk 'NR==2 {print $5}' | sed 's/%//'
}

# Function to get human readable disk info
get_disk_info() {
    df -h / | awk 'NR==2 {printf "Used: %s/%s (%s)", $3, $2, $5}'
}

# Function to run command or show what would run
run_or_show() {
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY RUN] Would run: $*"
    else
        echo "[CLEANUP] Running: $*"
        "$@"
    fi
}

echo "=== Docker Cleanup Script ==="
echo "Current disk usage: $(get_disk_info)"

CURRENT_USAGE=$(get_disk_usage)
echo "Disk usage: ${CURRENT_USAGE}%"

# Check if cleanup is needed
if [[ "$FORCE_CLEANUP" == "false" && "$CURRENT_USAGE" -lt "$DISK_THRESHOLD" ]]; then
    echo "✅ Disk usage (${CURRENT_USAGE}%) is below threshold (${DISK_THRESHOLD}%). No cleanup needed."
    exit 0
fi

echo "🧹 Starting Docker cleanup..."

# Stop and remove containers that are not running
echo "1. Cleaning up stopped containers..."
STOPPED_CONTAINERS=$(docker ps -aq --filter "status=exited" 2>/dev/null || true)
if [[ -n "$STOPPED_CONTAINERS" ]]; then
    run_or_show docker rm $STOPPED_CONTAINERS
else
    echo "   No stopped containers to remove"
fi

# Remove dangling images
echo "2. Removing dangling images..."
DANGLING_IMAGES=$(docker images -f "dangling=true" -q 2>/dev/null || true)
if [[ -n "$DANGLING_IMAGES" ]]; then
    run_or_show docker rmi $DANGLING_IMAGES
else
    echo "   No dangling images to remove"
fi

# Remove unused images (keep images from last 24 hours)
echo "3. Removing unused images older than 24 hours..."
run_or_show docker image prune -a -f --filter "until=24h"

# Remove unused volumes
echo "4. Removing unused volumes..."
run_or_show docker volume prune -f

# Remove unused networks
echo "5. Removing unused networks..."
run_or_show docker network prune -f

# Remove build cache (keep recent builds)
echo "6. Cleaning build cache..."
run_or_show docker builder prune -f --filter "until=24h"

# Final system prune for anything missed
echo "7. Final system cleanup..."
run_or_show docker system prune -f

echo ""
echo "=== Cleanup Complete ==="
echo "Disk usage after cleanup: $(get_disk_info)"

NEW_USAGE=$(get_disk_usage)
SAVED_SPACE=$((CURRENT_USAGE - NEW_USAGE))
echo "Space saved: ${SAVED_SPACE}% (${CURRENT_USAGE}% → ${NEW_USAGE}%)"

# Alert if still high usage
if [[ "$NEW_USAGE" -gt 90 ]]; then
    echo "⚠️  WARNING: Disk usage still high (${NEW_USAGE}%). Consider:"
    echo "   - Increasing server disk size"
    echo "   - Moving logs to external storage"
    echo "   - Implementing log rotation"
fi

echo "✅ Docker cleanup completed successfully"