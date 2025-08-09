#!/usr/bin/env bash
# Setup automated maintenance tasks
# Run this once to install cron jobs and monitoring

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="/opt/hackathon-server/sci-vis-ai-server"

echo "=== Setting up automated maintenance ==="

# Create maintenance directory
mkdir -p "$PROJECT_DIR/scripts/maintenance"

# Copy scripts to project directory
cp "$SCRIPT_DIR/docker-cleanup.sh" "$PROJECT_DIR/scripts/maintenance/"
cp "$SCRIPT_DIR/deploy-safe.sh" "$PROJECT_DIR/scripts/maintenance/"
chmod +x "$PROJECT_DIR/scripts/maintenance/"*.sh

# Create log directory
mkdir -p "$PROJECT_DIR/logs/maintenance"

# Create cron jobs
echo "Setting up cron jobs..."

# Create temporary cron file
TEMP_CRON=$(mktemp)

# Get existing cron jobs (if any)
crontab -l 2>/dev/null > "$TEMP_CRON" || true

# Add our maintenance jobs (if not already present)
if ! grep -q "docker-cleanup.sh" "$TEMP_CRON"; then
    echo "# Docker cleanup - runs daily at 2 AM" >> "$TEMP_CRON"
    echo "0 2 * * * $PROJECT_DIR/scripts/maintenance/docker-cleanup.sh --threshold=80 >> $PROJECT_DIR/logs/maintenance/cleanup.log 2>&1" >> "$TEMP_CRON"
fi

if ! grep -q "disk-monitor" "$TEMP_CRON"; then
    echo "# Disk space monitoring - runs every 4 hours" >> "$TEMP_CRON"
    echo "0 */4 * * * $PROJECT_DIR/scripts/maintenance/disk-monitor.sh >> $PROJECT_DIR/logs/maintenance/monitor.log 2>&1" >> "$TEMP_CRON"
fi

# Install the cron jobs
crontab "$TEMP_CRON"
rm "$TEMP_CRON"

echo "✅ Cron jobs installed:"
crontab -l | grep -E "(docker-cleanup|disk-monitor)"

# Create disk monitoring script
cat > "$PROJECT_DIR/scripts/maintenance/disk-monitor.sh" << 'EOF'
#!/usr/bin/env bash
# Disk space monitoring script

THRESHOLD=85
CRITICAL_THRESHOLD=95
PROJECT_DIR="/opt/hackathon-server/sci-vis-ai-server"

get_disk_usage() {
    df / | awk 'NR==2 {print $5}' | sed 's/%//'
}

get_disk_info() {
    df -h / | awk 'NR==2 {printf "%s/%s (%s)", $3, $2, $5}'
}

USAGE=$(get_disk_usage)
TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')

echo "[$TIMESTAMP] Disk usage: $(get_disk_info)"

if [[ "$USAGE" -gt "$CRITICAL_THRESHOLD" ]]; then
    echo "[$TIMESTAMP] CRITICAL: Disk usage ${USAGE}% > ${CRITICAL_THRESHOLD}%"
    echo "[$TIMESTAMP] Running emergency cleanup..."
    "$PROJECT_DIR/scripts/maintenance/docker-cleanup.sh" --force
elif [[ "$USAGE" -gt "$THRESHOLD" ]]; then
    echo "[$TIMESTAMP] WARNING: Disk usage ${USAGE}% > ${THRESHOLD}%"
    echo "[$TIMESTAMP] Running cleanup..."
    "$PROJECT_DIR/scripts/maintenance/docker-cleanup.sh" --threshold="$THRESHOLD"
else
    echo "[$TIMESTAMP] OK: Disk usage ${USAGE}% < ${THRESHOLD}%"
fi
EOF

chmod +x "$PROJECT_DIR/scripts/maintenance/disk-monitor.sh"

# Create log rotation config
cat > "/etc/logrotate.d/moleculens-maintenance" << EOF
$PROJECT_DIR/logs/maintenance/*.log {
    daily
    rotate 7
    compress
    delaycompress
    missingok
    notifempty
    create 644 root root
}
EOF

echo "✅ Log rotation configured"

# Create systemd service for monitoring (optional)
cat > "/etc/systemd/system/moleculens-monitor.service" << EOF
[Unit]
Description=Moleculens Disk Space Monitor
After=network.target

[Service]
Type=oneshot
ExecStart=$PROJECT_DIR/scripts/maintenance/disk-monitor.sh
User=root
StandardOutput=append:$PROJECT_DIR/logs/maintenance/monitor.log
StandardError=append:$PROJECT_DIR/logs/maintenance/monitor.log

[Install]
WantedBy=multi-user.target
EOF

# Create systemd timer
cat > "/etc/systemd/system/moleculens-monitor.timer" << EOF
[Unit]
Description=Run Moleculens Disk Monitor every 4 hours
Requires=moleculens-monitor.service

[Timer]
OnCalendar=*-*-* 00,04,08,12,16,20:00:00
Persistent=true

[Install]
WantedBy=timers.target
EOF

# Enable systemd timer (alternative to cron)
systemctl daemon-reload
systemctl enable moleculens-monitor.timer
systemctl start moleculens-monitor.timer

echo "✅ Systemd timer enabled as backup to cron"

# Create README for maintenance
cat > "$PROJECT_DIR/scripts/maintenance/README.md" << 'EOF'
# Maintenance Scripts

This directory contains automated maintenance scripts to prevent disk space issues.

## Scripts

- `docker-cleanup.sh` - Automated Docker cleanup with configurable thresholds
- `deploy-safe.sh` - Enhanced deployment script with disk space monitoring
- `disk-monitor.sh` - Continuous disk space monitoring
- `setup-maintenance.sh` - One-time setup script

## Automated Tasks

- **Daily cleanup**: Runs at 2 AM to clean Docker artifacts
- **Disk monitoring**: Checks every 4 hours and cleans if needed
- **Log rotation**: Keeps 7 days of maintenance logs

## Manual Usage

```bash
# Run cleanup manually
./docker-cleanup.sh --threshold=80

# Deploy safely with disk checks
./deploy-safe.sh

# Check current status
./disk-monitor.sh
```

## Monitoring

Check logs in `../../logs/maintenance/`:
- `cleanup.log` - Docker cleanup activities
- `monitor.log` - Disk space monitoring

## Configuration

Edit thresholds in scripts:
- `DISK_THRESHOLD=80` - Trigger cleanup at 80% usage
- `CRITICAL_THRESHOLD=95` - Emergency cleanup at 95%
EOF

echo ""
echo "=== Maintenance Setup Complete ==="
echo "✅ Scripts installed in: $PROJECT_DIR/scripts/maintenance/"
echo "✅ Cron jobs configured for automated cleanup"
echo "✅ Systemd timer enabled for monitoring"
echo "✅ Log rotation configured"
echo ""
echo "Next steps:"
echo "1. Test the scripts manually"
echo "2. Monitor logs in $PROJECT_DIR/logs/maintenance/"
echo "3. Adjust thresholds if needed"
echo ""
echo "The system will now automatically:"
echo "- Clean Docker artifacts daily at 2 AM"
echo "- Monitor disk space every 4 hours"
echo "- Prevent deployment failures due to disk space"