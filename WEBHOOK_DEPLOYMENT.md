# Webhook Deployment System

This document describes the automated deployment system that triggers when code is pushed to the main branch.

## 🚀 How It Works

When you push code to the `main` branch on GitHub:

1. **GitHub sends webhook** → Server receives push notification
2. **Webhook validates** → Verifies signature and checks branch
3. **Auto-deploy triggers** → Pulls latest code and rebuilds containers
4. **Health check** → Verifies deployment success
5. **Cleanup** → Removes old Docker artifacts
6. **Notification** → Logs deployment status

## 📋 System Components

### 1. Webhook Receiver (`scripts/webhook/webhook-receiver.py`)
- **Port**: 9000
- **Endpoint**: `/webhook`
- **Health Check**: `/health`
- **Security**: HMAC signature verification
- **Logging**: All events logged to `logs/webhook.log`

### 2. Auto-Deploy Script (`scripts/maintenance/auto-deploy.sh`)
- **Function**: Handles the actual deployment process
- **Features**: Disk space checks, rollback on failure, health verification
- **Logging**: Deployment logs to `logs/auto-deploy.log`
- **Timeout**: 10-minute maximum deployment time

### 3. Systemd Service (`moleculens-webhook.service`)
- **Auto-start**: Starts on boot
- **Auto-restart**: Restarts on failure
- **Resource limits**: 256MB RAM, 50% CPU
- **Security**: Restricted file system access

## 🔧 Setup Instructions

### Initial Setup (Run Once)

```bash
# SSH to server
ssh root@api.moleculens.com

# Navigate to project
cd /opt/hackathon-server/sci-vis-ai-server

# Run webhook setup
chmod +x scripts/webhook/setup-webhook.sh
./scripts/webhook/setup-webhook.sh
```

### GitHub Webhook Configuration

1. **Go to GitHub repository settings**:
   ```
   https://github.com/jakekinchen/moleculens-server/settings/hooks
   ```

2. **Click "Add webhook"**

3. **Configure webhook**:
   - **Payload URL**: `http://api.moleculens.com:9000/webhook`
   - **Content type**: `application/json`
   - **Secret**: `moleculens-webhook-secret-2024`
   - **Events**: Select "Just the push event"
   - **Active**: ✅ Checked

4. **Save webhook**

## 📊 Monitoring & Logs

### Service Status
```bash
# Check webhook service status
systemctl status moleculens-webhook.service

# Check if webhook is responding
curl http://localhost:9000/health

# View service logs
journalctl -u moleculens-webhook.service -f
```

### Deployment Logs
```bash
# View webhook logs
tail -f /opt/hackathon-server/sci-vis-ai-server/logs/webhook.log

# View deployment logs
tail -f /opt/hackathon-server/sci-vis-ai-server/logs/auto-deploy.log

# View recent deployments
grep "WEBHOOK TRIGGERED DEPLOYMENT" logs/auto-deploy.log | tail -5
```

### Real-time Monitoring
```bash
# Monitor all webhook activity
journalctl -u moleculens-webhook.service -f

# Monitor deployments in real-time
tail -f logs/auto-deploy.log

# Monitor both simultaneously
multitail logs/webhook.log logs/auto-deploy.log
```

## 🔍 Deployment Process Details

### Pre-deployment Checks
1. **Disk space verification** (fails if >95% usage)
2. **Docker service check**
3. **Git repository validation**
4. **Deployment lock check** (prevents concurrent deployments)

### Deployment Steps
1. **Backup current state** (for potential rollback)
2. **Pull latest changes** from origin/main
3. **Stop existing services** (graceful shutdown)
4. **Build new containers** (with timeout protection)
5. **Start new services**
6. **Health verification** (30-second timeout)
7. **Post-deployment cleanup**

### Rollback on Failure
- **Automatic rollback** to previous commit if deployment fails
- **Service restoration** attempts to restore previous working state
- **Notification system** alerts about rollback events

## 🛡️ Security Features

### Webhook Security
- **HMAC signature verification** using shared secret
- **IP validation** (GitHub webhook IPs only)
- **Path validation** (only `/webhook` endpoint accepts POST)
- **Payload size limits**

### System Security
- **Systemd security** (restricted file system access)
- **Resource limits** (prevents resource exhaustion)
- **User isolation** (runs as dedicated user)
- **Firewall rules** (only necessary ports open)

## 🚨 Troubleshooting

### Webhook Not Receiving Events

1. **Check service status**:
   ```bash
   systemctl status moleculens-webhook.service
   ```

2. **Check firewall**:
   ```bash
   ufw status
   netstat -tlnp | grep 9000
   ```

3. **Test webhook endpoint**:
   ```bash
   curl http://api.moleculens.com:9000/health
   ```

4. **Check GitHub webhook deliveries**:
   - Go to repository settings → Webhooks
   - Click on webhook → Recent Deliveries
   - Check for failed deliveries

### Deployment Failures

1. **Check deployment logs**:
   ```bash
   tail -50 logs/auto-deploy.log
   ```

2. **Check disk space**:
   ```bash
   df -h
   ```

3. **Check Docker status**:
   ```bash
   docker ps
   systemctl status docker
   ```

4. **Manual deployment**:
   ```bash
   ./scripts/maintenance/deploy-safe.sh
   ```

### Service Issues

1. **Restart webhook service**:
   ```bash
   systemctl restart moleculens-webhook.service
   ```

2. **Check service logs**:
   ```bash
   journalctl -u moleculens-webhook.service --since "1 hour ago"
   ```

3. **Validate configuration**:
   ```bash
   python3 -m py_compile scripts/webhook/webhook-receiver.py
   ```

## ⚙️ Configuration

### Environment Variables
```bash
# Webhook configuration
WEBHOOK_SECRET=moleculens-webhook-secret-2024
WEBHOOK_PORT=9000

# Deployment configuration
PROJECT_DIR=/opt/hackathon-server/sci-vis-ai-server
MAX_DEPLOY_TIME=600  # 10 minutes
```

### Customization Options

#### Change Webhook Port
```bash
# Edit systemd service
sudo systemctl edit moleculens-webhook.service

# Add override:
[Service]
Environment=WEBHOOK_PORT=8080

# Restart service
sudo systemctl restart moleculens-webhook.service
```

#### Change Webhook Secret
```bash
# Update systemd service
sudo systemctl edit moleculens-webhook.service

# Add override:
[Service]
Environment=WEBHOOK_SECRET=your-new-secret

# Update GitHub webhook configuration
# Restart service
sudo systemctl restart moleculens-webhook.service
```

#### Add Notifications
Edit `auto-deploy.sh` and uncomment/configure notification sections:
```bash
# Slack notifications
SLACK_WEBHOOK_URL="https://hooks.slack.com/services/..."

# Discord notifications  
DISCORD_WEBHOOK_URL="https://discord.com/api/webhooks/..."

# Email notifications
SMTP_SERVER="smtp.gmail.com"
```

## 📈 Performance & Scaling

### Resource Usage
- **Memory**: ~50MB for webhook receiver
- **CPU**: Minimal (only during deployments)
- **Disk**: Log files rotate automatically
- **Network**: Minimal (webhook events only)

### Scaling Considerations
- **Concurrent deployments**: Prevented by lock file
- **Multiple webhooks**: Service can handle multiple repositories
- **Load balancing**: Can run multiple webhook receivers on different ports

## 🔄 Maintenance

### Regular Tasks

#### Weekly
```bash
# Check service health
systemctl status moleculens-webhook.service

# Review deployment logs
tail -100 logs/auto-deploy.log | grep -E "(SUCCESS|FAILED)"

# Check disk usage
df -h
```

#### Monthly
```bash
# Rotate logs manually if needed
logrotate -f /etc/logrotate.d/moleculens-webhook

# Update webhook secret (optional)
# Review and update notification settings
```

### Log Rotation
Logs are automatically rotated by systemd journal, but you can configure custom rotation:

```bash
# Create logrotate config
cat > /etc/logrotate.d/moleculens-webhook << EOF
/opt/hackathon-server/sci-vis-ai-server/logs/*.log {
    daily
    rotate 14
    compress
    delaycompress
    missingok
    notifempty
    create 644 root root
}
EOF
```

## 🎯 Best Practices

### Development Workflow
1. **Test locally** before pushing to main
2. **Use feature branches** for experimental changes
3. **Monitor deployments** after pushing
4. **Check health** after automatic deployment

### Deployment Safety
1. **Small, frequent commits** are safer than large changes
2. **Database migrations** should be backward compatible
3. **Environment variables** should be managed separately
4. **Secrets** should never be in code

### Monitoring
1. **Set up alerts** for deployment failures
2. **Monitor resource usage** during deployments
3. **Keep deployment logs** for troubleshooting
4. **Test rollback procedures** periodically

---

## 🎉 Summary

The webhook deployment system provides:

- ✅ **Automatic deployment** on push to main
- ✅ **Security verification** with HMAC signatures
- ✅ **Rollback protection** on deployment failures
- ✅ **Health monitoring** and verification
- ✅ **Resource management** and cleanup
- ✅ **Comprehensive logging** for troubleshooting
- ✅ **Service reliability** with auto-restart

**Your deployments are now fully automated and safe!** 🚀

Push to main → Automatic deployment → Service running with latest code!