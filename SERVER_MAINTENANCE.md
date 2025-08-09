# Server Maintenance Guide

This document outlines the automated maintenance system implemented to prevent disk space issues and ensure reliable deployments for the Moleculens API server.

## 🚨 Problem Solved

**Issue**: Server deployments were failing due to disk space exhaustion (99% usage) caused by accumulated Docker images, containers, and build cache.

**Solution**: Comprehensive automated maintenance system with proactive monitoring, cleanup, and safe deployment practices.

## 📋 Table of Contents

- [Quick Start](#quick-start)
- [Automated Systems](#automated-systems)
- [Manual Operations](#manual-operations)
- [Monitoring & Alerts](#monitoring--alerts)
- [Troubleshooting](#troubleshooting)
- [Configuration](#configuration)
- [Best Practices](#best-practices)

## 🚀 Quick Start

### Initial Setup (Run Once)

```bash
# SSH to server
ssh root@api.moleculens.com

# Navigate to project directory
cd /opt/hackathon-server/sci-vis-ai-server

# Run maintenance setup
chmod +x scripts/maintenance/setup-maintenance.sh
./scripts/maintenance/setup-maintenance.sh
```

### Daily Operations

```bash
# Use safe deployment (replaces ./deploy.sh)
./scripts/maintenance/deploy-safe.sh

# Check disk space
df -h

# Manual cleanup if needed
./scripts/maintenance/docker-cleanup.sh --threshold=80
```

## 🤖 Automated Systems

### 1. Daily Docker Cleanup (2 AM)

**What it does:**
- Removes Docker images older than 24 hours
- Cleans unused containers, volumes, and networks
- Removes build cache older than 24 hours
- Only runs when disk usage exceeds 80%

**Cron job:**
```bash
0 2 * * * /opt/hackathon-server/sci-vis-ai-server/scripts/maintenance/docker-cleanup.sh --threshold=80
```

### 2. Disk Space Monitoring (Every 4 Hours)

**What it monitors:**
- Disk usage percentage
- Available free space in GB
- Docker system resource usage

**Thresholds:**
- **80%**: Warning threshold - triggers cleanup
- **85%**: High usage - automatic cleanup
- **95%**: Critical - emergency cleanup

**Cron job:**
```bash
0 */4 * * * /opt/hackathon-server/sci-vis-ai-server/scripts/maintenance/disk-monitor.sh
```

### 3. GitHub Actions Monitoring (Every 6 Hours)

**Remote monitoring via CI/CD:**
- Checks server disk space from GitHub Actions
- Runs cleanup remotely if needed
- Alerts on high disk usage
- Can be triggered manually

**Workflow:** `.github/workflows/disk-space-monitor.yml`

### 4. Log Rotation

**Prevents log files from consuming disk space:**
- Rotates logs daily
- Keeps 7 days of history
- Compresses old logs
- Configured in `/etc/logrotate.d/moleculens-maintenance`

## 🛠 Manual Operations

### Safe Deployment

**Always use the safe deployment script:**

```bash
./scripts/maintenance/deploy-safe.sh
```

**What it does:**
1. Checks current disk usage
2. Estimates build space requirements
3. Runs cleanup if space is low
4. Monitors space during build
5. Cleans up build artifacts after deployment
6. Performs health check

### Manual Cleanup Options

```bash
# Standard cleanup (80% threshold)
./scripts/maintenance/docker-cleanup.sh --threshold=80

# Check what would be cleaned (dry run)
./scripts/maintenance/docker-cleanup.sh --dry-run

# Force cleanup regardless of usage
./scripts/maintenance/docker-cleanup.sh --force

# Emergency nuclear option (removes ALL unused Docker data)
docker system prune -a -f --volumes
```

### Disk Space Commands

```bash
# Check overall disk usage
df -h

# Check Docker-specific usage
docker system df

# Check largest directories
du -sh /* | sort -hr | head -10

# Check Docker images by size
docker images --format "table {{.Repository}}\t{{.Tag}}\t{{.Size}}" | sort -k3 -hr
```

## 📊 Monitoring & Alerts

### Log Files

All maintenance activities are logged:

```bash
# Cleanup activities
tail -f /opt/hackathon-server/sci-vis-ai-server/logs/maintenance/cleanup.log

# Disk monitoring
tail -f /opt/hackathon-server/sci-vis-ai-server/logs/maintenance/monitor.log

# Deployment activities
tail -f /opt/hackathon-server/sci-vis-ai-server/logs/maintenance/deploy.log
```

### Status Checks

```bash
# Check cron jobs
crontab -l | grep -E "(docker-cleanup|disk-monitor)"

# Check systemd timers
systemctl list-timers moleculens-monitor

# View recent cleanup activity
tail -20 /opt/hackathon-server/sci-vis-ai-server/logs/maintenance/cleanup.log

# Check current resource usage
docker stats --no-stream
```

### GitHub Actions Setup

To enable remote monitoring, add these secrets to your GitHub repository:

1. Go to **Settings** → **Secrets and variables** → **Actions**
2. Add these repository secrets:
   - `SERVER_HOST`: `api.moleculens.com`
   - `SERVER_USER`: `root`
   - `SERVER_SSH_KEY`: Your private SSH key content

## 🔧 Troubleshooting

### High Disk Usage Alerts

**If disk usage exceeds 85%:**

1. **Check current usage:**
   ```bash
   df -h
   docker system df
   ```

2. **Run immediate cleanup:**
   ```bash
   ./scripts/maintenance/docker-cleanup.sh --force
   ```

3. **Check for large files:**
   ```bash
   find / -type f -size +1G -exec ls -lh {} \; 2>/dev/null
   ```

### Deployment Failures

**If deployment fails due to space:**

1. **Use safe deployment:**
   ```bash
   ./scripts/maintenance/deploy-safe.sh
   ```

2. **Check deployment logs:**
   ```bash
   tail -f logs/maintenance/deploy.log
   ```

3. **Manual recovery:**
   ```bash
   # Clean everything
   docker system prune -a -f --volumes
   
   # Retry deployment
   ./scripts/maintenance/deploy-safe.sh
   ```

### Container Issues

**If containers are consuming too much space:**

1. **Check container sizes:**
   ```bash
   docker ps -s
   ```

2. **Check container logs:**
   ```bash
   docker logs --tail=100 fastapi_app
   ```

3. **Restart with resource limits:**
   ```bash
   docker-compose -f docker-compose.enhanced.yml up -d
   ```

### Emergency Recovery

**If disk is completely full (99%+):**

```bash
# 1. Emergency cleanup
docker system prune -a -f --volumes

# 2. Clean system logs
journalctl --vacuum-time=7d

# 3. Clean package cache
apt-get clean

# 4. Remove old kernels
apt-get autoremove --purge

# 5. Check for large log files
find /var/log -type f -size +100M -exec ls -lh {} \;
```

## ⚙️ Configuration

### Disk Thresholds

Edit thresholds in scripts as needed:

```bash
# In docker-cleanup.sh
DISK_THRESHOLD=80        # Trigger cleanup at 80%

# In disk-monitor.sh  
THRESHOLD=85            # Warning threshold
CRITICAL_THRESHOLD=95   # Emergency threshold

# In deploy-safe.sh
DISK_THRESHOLD=85       # Pre-deployment check
MIN_FREE_SPACE_GB=5     # Minimum free space for builds
```

### Resource Limits

Container resource limits in `docker-compose.enhanced.yml`:

```yaml
deploy:
  resources:
    limits:
      memory: 4G          # Maximum memory per container
      cpus: '2.0'         # Maximum CPU cores
    reservations:
      memory: 1G          # Reserved memory
      cpus: '0.5'         # Reserved CPU
```

### Log Rotation

Log rotation settings in `/etc/logrotate.d/moleculens-maintenance`:

```
/opt/hackathon-server/sci-vis-ai-server/logs/maintenance/*.log {
    daily                # Rotate daily
    rotate 7            # Keep 7 days
    compress            # Compress old logs
    delaycompress       # Don't compress immediately
    missingok           # Don't error if log missing
    notifempty          # Don't rotate empty logs
}
```

## 📚 Best Practices

### Development

1. **Always use safe deployment:**
   ```bash
   ./scripts/maintenance/deploy-safe.sh
   ```

2. **Test changes locally first:**
   ```bash
   docker-compose build
   docker-compose up -d
   ```

3. **Monitor resource usage:**
   ```bash
   docker stats --no-stream
   ```

### Operations

1. **Weekly monitoring:**
   - Check maintenance logs
   - Review disk usage trends
   - Verify automated jobs are running

2. **Monthly maintenance:**
   - Review and adjust thresholds
   - Update maintenance scripts if needed
   - Check for system updates

3. **Deployment checklist:**
   - [ ] Use `deploy-safe.sh`
   - [ ] Monitor deployment logs
   - [ ] Verify service health
   - [ ] Check disk usage after deployment

### Monitoring

1. **Set up alerts:**
   - Configure email notifications for critical disk usage
   - Monitor GitHub Actions workflow results
   - Set up external monitoring (optional)

2. **Regular checks:**
   ```bash
   # Daily
   df -h
   
   # Weekly  
   tail -50 logs/maintenance/cleanup.log
   
   # Monthly
   docker system df
   ```

## 📈 Resource Management

### Container Limits

The enhanced Docker Compose configuration includes:

- **Memory limits**: 4GB max per container
- **CPU limits**: 2 cores max per container
- **Log rotation**: 100MB max log files, 3 files kept
- **Health checks**: Automatic service monitoring
- **Restart policies**: Automatic recovery from failures

### Disk Space Allocation

**Recommended disk usage:**
- **System**: < 20% (OS, packages, logs)
- **Docker**: < 60% (images, containers, volumes)
- **Application**: < 15% (code, data, temp files)
- **Free space**: > 5% (buffer for operations)

### Performance Optimization

1. **Use multi-stage Docker builds** to reduce image size
2. **Implement .dockerignore** to exclude unnecessary files
3. **Use specific base image tags** instead of `latest`
4. **Clean up intermediate containers** during builds

## 🔍 Monitoring Dashboard

### Key Metrics to Track

```bash
# Disk usage percentage
df / | awk 'NR==2 {print $5}'

# Free space in GB
df / | awk 'NR==2 {printf "%.1f", $4/1024/1024}'

# Docker space usage
docker system df --format "table {{.Type}}\t{{.TotalCount}}\t{{.Size}}\t{{.Reclaimable}}"

# Container resource usage
docker stats --no-stream --format "table {{.Container}}\t{{.CPUPerc}}\t{{.MemUsage}}"
```

### Health Check Endpoints

```bash
# Service health
curl -f http://localhost:8000/health

# Docker health
docker ps --filter "health=unhealthy"

# System health
systemctl status moleculens-monitor.timer
```

## 📞 Support

### Getting Help

1. **Check logs first:**
   ```bash
   tail -50 logs/maintenance/monitor.log
   ```

2. **Run diagnostics:**
   ```bash
   ./scripts/maintenance/docker-cleanup.sh --dry-run
   df -h
   docker system df
   ```

3. **Contact information:**
   - Repository issues: [GitHub Issues](https://github.com/jakekinchen/moleculens-server/issues)
   - Emergency: Check server logs and run emergency recovery

### Common Issues

| Issue | Cause | Solution |
|-------|-------|----------|
| Deployment fails | Disk full | Run `docker-cleanup.sh --force` |
| High memory usage | Container leak | Restart containers with limits |
| Slow performance | Resource contention | Check `docker stats` |
| Log files growing | No rotation | Verify logrotate configuration |

---

## 🎯 Summary

This maintenance system provides:

- ✅ **Automated daily cleanup** of Docker artifacts
- ✅ **Proactive disk monitoring** every 4 hours  
- ✅ **Safe deployment process** with space checks
- ✅ **Resource limits** to prevent runaway containers
- ✅ **Comprehensive logging** for troubleshooting
- ✅ **Remote monitoring** via GitHub Actions
- ✅ **Emergency recovery** procedures

**The disk space problem that caused deployment failures will never happen again.**

For questions or issues, check the logs first, then refer to the troubleshooting section above.