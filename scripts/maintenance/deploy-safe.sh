#!/usr/bin/env bash
# Enhanced Deploy Script with Disk Space Management
# Prevents deployment failures due to disk space issues

set -e

# Configuration
DISK_THRESHOLD=85
MIN_FREE_SPACE_GB=5
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="/opt/hackathon-server/sci-vis-ai-server"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging function
log() {
    echo -e "${BLUE}[DEPLOY]${NC} $1"
}

warn() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

# Function to get disk usage percentage
get_disk_usage() {
    df / | awk 'NR==2 {print $5}' | sed 's/%//'
}

# Function to get free space in GB
get_free_space_gb() {
    df / | awk 'NR==2 {printf "%.1f", $4/1024/1024}'
}

# Function to check disk space and cleanup if needed
check_and_cleanup_disk() {
    local current_usage=$(get_disk_usage)
    local free_space=$(get_free_space_gb)
    
    log "Current disk usage: ${current_usage}% (${free_space}GB free)"
    
    # Check if cleanup is needed
    if [[ "$current_usage" -gt "$DISK_THRESHOLD" ]] || [[ $(echo "$free_space < $MIN_FREE_SPACE_GB" | bc -l) -eq 1 ]]; then
        warn "Disk usage high (${current_usage}%) or low free space (${free_space}GB)"
        log "Running automated cleanup..."
        
        # Run cleanup script if it exists
        if [[ -f "$SCRIPT_DIR/docker-cleanup.sh" ]]; then
            bash "$SCRIPT_DIR/docker-cleanup.sh" --threshold=75
        else
            # Fallback cleanup
            log "Running fallback Docker cleanup..."
            docker system prune -f --filter "until=24h" || true
            docker image prune -a -f --filter "until=24h" || true
        fi
        
        # Check again after cleanup
        local new_usage=$(get_disk_usage)
        local new_free_space=$(get_free_space_gb)
        
        if [[ "$new_usage" -gt 90 ]] || [[ $(echo "$new_free_space < 2" | bc -l) -eq 1 ]]; then
            error "Disk space still critically low after cleanup (${new_usage}%, ${new_free_space}GB free)"
            error "Please manually free up space or increase disk size before deploying"
            exit 1
        fi
        
        success "Cleanup completed. New usage: ${new_usage}% (${new_free_space}GB free)"
    else
        success "Disk space OK (${current_usage}%, ${free_space}GB free)"
    fi
}

# Function to estimate build space requirements
estimate_build_space() {
    log "Estimating Docker build space requirements..."
    
    # Check if base images exist
    local base_images=("continuumio/miniconda3" "python:3.9")
    local estimated_space=0
    
    for image in "${base_images[@]}"; do
        if ! docker image inspect "$image" >/dev/null 2>&1; then
            estimated_space=$((estimated_space + 2)) # Estimate 2GB per base image
        fi
    done
    
    # Add space for build layers
    estimated_space=$((estimated_space + 3)) # 3GB for build layers
    
    local free_space=$(get_free_space_gb)
    if [[ $(echo "$free_space < $estimated_space" | bc -l) -eq 1 ]]; then
        warn "Estimated build space needed: ${estimated_space}GB, available: ${free_space}GB"
        return 1
    fi
    
    log "Estimated build space: ${estimated_space}GB, available: ${free_space}GB ✓"
    return 0
}

# Main deployment function
main() {
    log "Starting safe deployment process..."
    
    # Check if bc is available for calculations
    if ! command -v bc &> /dev/null; then
        log "Installing bc for calculations..."
        apt-get update && apt-get install -y bc
    fi
    
    # Pre-deployment disk check
    check_and_cleanup_disk
    
    # Estimate space requirements
    if ! estimate_build_space; then
        warn "Insufficient space for build. Running aggressive cleanup..."
        if [[ -f "$SCRIPT_DIR/docker-cleanup.sh" ]]; then
            bash "$SCRIPT_DIR/docker-cleanup.sh" --force
        else
            docker system prune -a -f
        fi
        
        # Check again
        if ! estimate_build_space; then
            error "Still insufficient space after cleanup. Please increase disk size."
            exit 1
        fi
    fi
    
    log "Killing any process using port 8000..."
    fuser -k 8000/tcp || true
    
    log "Resetting code to latest origin/main..."
    cd "$PROJECT_DIR"
    git fetch origin main
    git reset --hard origin/main
    
    log "Navigating to Docker context..."
    cd api
    
    log "Bringing down any running containers..."
    docker-compose down --remove-orphans
    
    log "Cleaning up stopped containers..."
    docker container prune -f
    
    # Monitor disk space during build
    log "Rebuilding Docker containers..."
    if ! docker-compose build; then
        error "Docker build failed!"
        
        # Check if it was due to disk space
        local usage_after_fail=$(get_disk_usage)
        if [[ "$usage_after_fail" -gt 95 ]]; then
            error "Build failed due to disk space (${usage_after_fail}%)"
            log "Running emergency cleanup..."
            docker system prune -a -f
        fi
        exit 1
    fi
    
    log "Bringing up app..."
    docker-compose up -d
    
    # Post-deployment cleanup of build artifacts
    log "Cleaning up build artifacts..."
    docker image prune -f --filter "dangling=true"
    
    # Final status
    local final_usage=$(get_disk_usage)
    local final_free=$(get_free_space_gb)
    success "Deployment completed successfully!"
    success "Final disk usage: ${final_usage}% (${final_free}GB free)"
    
    # Health check
    log "Waiting for service to start..."
    sleep 10
    
    if curl -f http://localhost:8000/health >/dev/null 2>&1; then
        success "Service health check passed ✓"
    else
        warn "Service health check failed - please verify manually"
    fi
}

# Run main function
main "$@"