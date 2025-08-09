#!/usr/bin/env bash
# Automatic Deployment Script for Webhook Triggers
# This script is called by the webhook receiver when code is pushed to main

set -e

# Configuration
PROJECT_DIR="/opt/hackathon-server/sci-vis-ai-server"
LOG_FILE="$PROJECT_DIR/logs/auto-deploy.log"
LOCK_FILE="/tmp/moleculens-deploy.lock"
MAX_DEPLOY_TIME=600  # 10 minutes timeout

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging function
log() {
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo -e "${BLUE}[$timestamp]${NC} $1" | tee -a "$LOG_FILE"
}

warn() {
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo -e "${YELLOW}[$timestamp] WARNING:${NC} $1" | tee -a "$LOG_FILE"
}

error() {
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo -e "${RED}[$timestamp] ERROR:${NC} $1" | tee -a "$LOG_FILE"
}

success() {
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo -e "${GREEN}[$timestamp] SUCCESS:${NC} $1" | tee -a "$LOG_FILE"
}

# Function to cleanup on exit
cleanup() {
    if [[ -f "$LOCK_FILE" ]]; then
        rm -f "$LOCK_FILE"
    fi
}

# Function to send notification (placeholder for future integrations)
send_notification() {
    local status="$1"
    local message="$2"
    
    # Log the notification
    log "NOTIFICATION [$status]: $message"
    
    # Future: Add Slack, Discord, email notifications here
    # curl -X POST -H 'Content-type: application/json' \
    #   --data "{\"text\":\"🚀 Moleculens Deploy [$status]: $message\"}" \
    #   "$SLACK_WEBHOOK_URL"
}

# Function to check if deployment is already running
check_deployment_lock() {
    if [[ -f "$LOCK_FILE" ]]; then
        local lock_pid=$(cat "$LOCK_FILE" 2>/dev/null || echo "")
        if [[ -n "$lock_pid" ]] && kill -0 "$lock_pid" 2>/dev/null; then
            error "Deployment already in progress (PID: $lock_pid)"
            exit 1
        else
            warn "Stale lock file found, removing..."
            rm -f "$LOCK_FILE"
        fi
    fi
    
    # Create lock file
    echo $$ > "$LOCK_FILE"
}

# Function to get disk usage
get_disk_usage() {
    df / | awk 'NR==2 {print $5}' | sed 's/%//'
}

# Function to get free space in GB
get_free_space_gb() {
    df / | awk 'NR==2 {printf "%.1f", $4/1024/1024}'
}

# Function to check prerequisites
check_prerequisites() {
    log "Checking deployment prerequisites..."
    
    # Check if we're in the right directory
    if [[ ! -d "$PROJECT_DIR/.git" ]]; then
        error "Not in a git repository: $PROJECT_DIR"
        exit 1
    fi
    
    # Check disk space
    local disk_usage=$(get_disk_usage)
    local free_space=$(get_free_space_gb)
    
    log "Disk usage: ${disk_usage}% (${free_space}GB free)"
    
    if [[ "$disk_usage" -gt 90 ]]; then
        warn "High disk usage detected (${disk_usage}%)"
        log "Running emergency cleanup..."
        
        if [[ -f "$PROJECT_DIR/scripts/maintenance/docker-cleanup.sh" ]]; then
            "$PROJECT_DIR/scripts/maintenance/docker-cleanup.sh" --force
        else
            docker system prune -a -f || true
        fi
        
        # Check again after cleanup
        disk_usage=$(get_disk_usage)
        if [[ "$disk_usage" -gt 95 ]]; then
            error "Insufficient disk space after cleanup (${disk_usage}%)"
            send_notification "FAILED" "Deployment failed: insufficient disk space (${disk_usage}%)"
            exit 1
        fi
    fi
    
    # Check if Docker is running
    if ! docker info >/dev/null 2>&1; then
        error "Docker is not running"
        exit 1
    fi
    
    success "Prerequisites check passed"
}

# Function to backup current state
backup_current_state() {
    log "Creating backup of current state..."
    
    # Get current commit hash
    local current_commit=$(git rev-parse HEAD 2>/dev/null || echo "unknown")
    
    # Store in environment for potential rollback
    export BACKUP_COMMIT="$current_commit"
    
    log "Current commit: $current_commit"
}

# Function to pull latest changes
pull_latest_changes() {
    log "Pulling latest changes from origin/main..."
    
    cd "$PROJECT_DIR"
    
    # Fetch latest changes
    if ! git fetch origin main; then
        error "Failed to fetch from origin"
        exit 1
    fi
    
    # Get the latest commit info
    local latest_commit=$(git rev-parse origin/main)
    local current_commit=$(git rev-parse HEAD)
    
    if [[ "$latest_commit" == "$current_commit" ]]; then
        log "Already up to date with origin/main"
        return 0
    fi
    
    # Reset to latest
    if ! git reset --hard origin/main; then
        error "Failed to reset to origin/main"
        exit 1
    fi
    
    # Log the changes
    local new_commit=$(git rev-parse HEAD)
    local commit_message=$(git log -1 --pretty=format:"%s" "$new_commit")
    local author=$(git log -1 --pretty=format:"%an" "$new_commit")
    
    log "Updated to commit: ${new_commit:0:8}"
    log "Author: $author"
    log "Message: $commit_message"
    
    success "Successfully pulled latest changes"
}

# Function to deploy application
deploy_application() {
    log "Starting application deployment..."
    
    cd "$PROJECT_DIR"
    
    # Kill any processes using port 8000
    log "Stopping services on port 8000..."
    fuser -k 8000/tcp || true
    sleep 2
    
    # Navigate to API directory
    cd api
    
    # Stop existing containers
    log "Stopping existing containers..."
    cd /opt/hackathon-server/sci-vis-ai-server/api && docker-compose down --remove-orphans || true
    
    # Clean up stopped containers
    log "Cleaning up stopped containers..."
    docker container prune -f || true
    
    # Build new containers
    log "Building new containers..."
    if ! timeout $MAX_DEPLOY_TIME cd /opt/hackathon-server/sci-vis-ai-server/api && docker-compose build; then
        error "Docker build failed or timed out"
        return 1
    fi
    
    # Start new containers
    log "Starting new containers..."
    if ! cd /opt/hackathon-server/sci-vis-ai-server/api && docker-compose up -d; then
        error "Failed to start containers"
        return 1
    fi
    
    # Wait for service to be ready
    log "Waiting for service to be ready..."
    local max_attempts=30
    local attempt=0
    
    while [[ $attempt -lt $max_attempts ]]; do
        if curl -f http://localhost:8000/health >/dev/null 2>&1; then
            success "Service is healthy and ready"
            return 0
        fi
        
        attempt=$((attempt + 1))
        log "Waiting for service... (attempt $attempt/$max_attempts)"
        sleep 2
    done
    
    error "Service failed to become healthy within timeout"
    return 1
}

# Function to verify deployment
verify_deployment() {
    log "Verifying deployment..."
    
    # Check if containers are running
    local running_containers=$(docker ps --filter "name=fastapi_app" --format "{{.Names}}" | wc -l)
    if [[ "$running_containers" -eq 0 ]]; then
        error "No containers are running"
        return 1
    fi
    
    # Check service health
    if ! curl -f http://localhost:8000/health >/dev/null 2>&1; then
        error "Health check failed"
        return 1
    fi
    
    # Check if we can access the API docs
    if ! curl -f http://localhost:8000/docs >/dev/null 2>&1; then
        warn "API docs endpoint not accessible"
    fi
    
    # Log container status
    log "Container status:"
    docker ps --filter "name=fastapi_app" --format "table {{.Names}}\t{{.Status}}\t{{.Ports}}" | tee -a "$LOG_FILE"
    
    success "Deployment verification passed"
    return 0
}

# Function to cleanup after deployment
post_deploy_cleanup() {
    log "Running post-deployment cleanup..."
    
    # Clean up dangling images
    docker image prune -f >/dev/null 2>&1 || true
    
    # Clean up old build cache
    docker builder prune -f --filter "until=24h" >/dev/null 2>&1 || true
    
    # Log final disk usage
    local final_usage=$(get_disk_usage)
    local final_free=$(get_free_space_gb)
    log "Final disk usage: ${final_usage}% (${final_free}GB free)"
    
    success "Post-deployment cleanup completed"
}

# Function to rollback on failure
rollback_deployment() {
    error "Deployment failed, attempting rollback..."
    
    if [[ -n "$BACKUP_COMMIT" && "$BACKUP_COMMIT" != "unknown" ]]; then
        cd "$PROJECT_DIR"
        
        log "Rolling back to commit: $BACKUP_COMMIT"
        if git reset --hard "$BACKUP_COMMIT"; then
            log "Redeploying previous version..."
            cd api
            cd /opt/hackathon-server/sci-vis-ai-server/api && docker-compose down --remove-orphans || true
            cd /opt/hackathon-server/sci-vis-ai-server/api && docker-compose up -d || true
            
            # Wait a bit and check
            sleep 10
            if curl -f http://localhost:8000/health >/dev/null 2>&1; then
                warn "Rollback successful - service is running previous version"
                send_notification "ROLLBACK" "Deployment failed, rolled back to previous version"
            else
                error "Rollback failed - manual intervention required"
                send_notification "CRITICAL" "Deployment and rollback both failed - manual intervention required"
            fi
        else
            error "Failed to rollback to previous commit"
            send_notification "CRITICAL" "Rollback failed - manual intervention required"
        fi
    else
        error "No backup commit available for rollback"
        send_notification "CRITICAL" "Deployment failed and no rollback available"
    fi
}

# Main deployment function
main() {
    # Set up logging
    mkdir -p "$(dirname "$LOG_FILE")"
    
    # Set trap for cleanup
    trap cleanup EXIT
    
    # Log webhook information if available
    if [[ -n "$WEBHOOK_COMMIT_HASH" ]]; then
        log "=== WEBHOOK TRIGGERED DEPLOYMENT ==="
        log "Commit: $WEBHOOK_COMMIT_HASH"
        log "Author: $WEBHOOK_AUTHOR"
        log "Message: $WEBHOOK_COMMIT_MESSAGE"
        log "Timestamp: $WEBHOOK_TIMESTAMP"
        log "======================================="
    else
        log "=== MANUAL DEPLOYMENT ==="
    fi
    
    # Check if another deployment is running
    check_deployment_lock
    
    # Start deployment process
    local start_time=$(date +%s)
    
    if check_prerequisites && \
       backup_current_state && \
       pull_latest_changes && \
       deploy_application && \
       verify_deployment; then
        
        # Success path
        post_deploy_cleanup
        
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        
        success "Deployment completed successfully in ${duration} seconds"
        
        # Get final commit info
        local final_commit=$(git rev-parse HEAD 2>/dev/null || echo "unknown")
        local commit_message=$(git log -1 --pretty=format:"%s" "$final_commit" 2>/dev/null || echo "unknown")
        
        send_notification "SUCCESS" "Deployment completed: ${final_commit:0:8} - $commit_message"
        
    else
        # Failure path
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        
        error "Deployment failed after ${duration} seconds"
        rollback_deployment
        exit 1
    fi
}

# Run main function
main "$@"
