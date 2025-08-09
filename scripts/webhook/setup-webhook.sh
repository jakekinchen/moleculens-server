#!/usr/bin/env bash
# Setup script for GitHub webhook system

set -e

PROJECT_DIR="/opt/hackathon-server/sci-vis-ai-server"
WEBHOOK_SECRET="moleculens-webhook-secret-2025"

echo "🔗 Setting up GitHub webhook system..."

# Create webhook directory
mkdir -p "$PROJECT_DIR/scripts/webhook"
mkdir -p "$PROJECT_DIR/logs"

# Make scripts executable
chmod +x "$PROJECT_DIR/scripts/webhook/webhook-receiver.py"
chmod +x "$PROJECT_DIR/scripts/maintenance/auto-deploy.sh"

# Create systemd service
cat > /etc/systemd/system/moleculens-webhook.service << 'SYSTEMD_EOF'
[Unit]
Description=Moleculens GitHub Webhook Receiver
After=network.target
Wants=network.target

[Service]
Type=simple
User=root
Group=root
WorkingDirectory=/opt/hackathon-server/sci-vis-ai-server
ExecStart=/usr/bin/python3 /opt/hackathon-server/sci-vis-ai-server/scripts/webhook/webhook-receiver.py
Restart=always
RestartSec=10
StandardOutput=journal
StandardError=journal

# Environment variables
Environment=WEBHOOK_SECRET=moleculens-webhook-secret-2024
Environment=WEBHOOK_PORT=9000
Environment=PYTHONUNBUFFERED=1

# Security settings
NoNewPrivileges=true
PrivateTmp=true
ProtectSystem=strict
ReadWritePaths=/opt/hackathon-server/sci-vis-ai-server
ReadWritePaths=/tmp
ReadWritePaths=/var/run/docker.sock

# Resource limits
MemoryMax=256M
CPUQuota=50%

[Install]
WantedBy=multi-user.target
SYSTEMD_EOF

# Reload systemd and enable service
systemctl daemon-reload
systemctl enable moleculens-webhook.service
systemctl start moleculens-webhook.service

# Check if service is running
sleep 2
if systemctl is-active --quiet moleculens-webhook.service; then
    echo "✅ Webhook service is running"
else
    echo "❌ Webhook service failed to start"
    systemctl status moleculens-webhook.service
    exit 1
fi

# Test webhook endpoint
echo "🧪 Testing webhook endpoint..."
if curl -f http://localhost:9000/health >/dev/null 2>&1; then
    echo "✅ Webhook endpoint is responding"
else
    echo "❌ Webhook endpoint is not responding"
    exit 1
fi

# Open firewall port if ufw is active
if command -v ufw >/dev/null && ufw status | grep -q "Status: active"; then
    echo "🔥 Opening firewall port 9000..."
    ufw allow 9000/tcp
fi

echo ""
echo "✅ Webhook system setup complete!"
echo ""
echo "📋 Next steps:"
echo "1. Add webhook to GitHub repository:"
echo "   - Go to: https://github.com/jakekinchen/moleculens-server/settings/hooks"
echo "   - Click 'Add webhook'"
echo "   - Payload URL: http://api.moleculens.com:9000/webhook"
echo "   - Content type: application/json"
echo "   - Secret: $WEBHOOK_SECRET"
echo "   - Events: Just the push event"
echo ""
echo "2. Test the webhook:"
echo "   - Make a commit and push to main branch"
echo "   - Check logs: journalctl -u moleculens-webhook.service -f"
echo "   - Check deploy logs: tail -f $PROJECT_DIR/logs/auto-deploy.log"
echo ""
echo "🔍 Monitoring commands:"
echo "   - Service status: systemctl status moleculens-webhook.service"
echo "   - Service logs: journalctl -u moleculens-webhook.service -f"
echo "   - Webhook health: curl http://localhost:9000/health"
echo "   - Deploy logs: tail -f $PROJECT_DIR/logs/auto-deploy.log"
