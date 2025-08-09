# Server Operations Guide

This document provides comprehensive instructions for all Git/GitHub operations and server management tasks for the Moleculens API server.

## 📋 Table of Contents

- [Server Information](#server-information)
- [Git & GitHub Operations](#git--github-operations)
- [Deployment Operations](#deployment-operations)
- [Docker Operations](#docker-operations)
- [Service Management](#service-management)
- [Monitoring & Debugging](#monitoring--debugging)
- [File System Operations](#file-system-operations)
- [Security & Access](#security--access)
- [Emergency Procedures](#emergency-procedures)
- [Best Practices](#best-practices)

## 🖥️ Server Information

### Server Details
- **Host**: `api.moleculens.com`
- **IP**: `165.232.151.162`
- **OS**: Ubuntu 24.04 LTS
- **Architecture**: x86_64
- **Resources**: 2 vCPU, 8GB RAM, 28GB SSD

### Key Directories
```bash
/opt/hackathon-server/sci-vis-ai-server/    # Main project directory
/opt/hackathon-server/sci-vis-ai-server/api/    # API source code
/opt/hackathon-server/sci-vis-ai-server/logs/   # Application logs
/var/log/                                   # System logs
/etc/ssh/                                   # SSH configuration
```

### Services Running
- **FastAPI Application**: Port 8000
- **Docker**: Container management
- **Nginx**: Reverse proxy (if configured)
- **SSH**: Port 22
- **Cron**: Automated maintenance tasks

## 🔧 Git & GitHub Operations

### SSH Access to Server
```bash
# Connect to server
ssh root@api.moleculens.com

# Or with specific key
ssh -i ~/.ssh/your-key root@api.moleculens.com
```

### Repository Management

#### Basic Git Operations
```bash
# Navigate to project
cd /opt/hackathon-server/sci-vis-ai-server

# Check current status
git status

# View recent commits
git log --oneline -10

# Check current branch
git branch

# View remote information
git remote -v
```

#### Pulling Latest Changes
```bash
# Fetch latest changes
git fetch origin main

# Pull and merge (safe)
git pull origin main

# Hard reset to latest (destructive - loses local changes)
git reset --hard origin/main
```

#### Making Changes
```bash
# Check what files changed
git status
git diff

# Add specific files
git add filename.py
git add directory/

# Add all changes
git add .

# Commit with message
git commit -m "feat: add new feature description"

# Push to GitHub
git push origin main
```

#### Branch Operations
```bash
# Create new branch
git checkout -b feature/new-feature

# Switch branches
git checkout main
git checkout feature/new-feature

# List all branches
git branch -a

# Delete branch
git branch -d feature/old-feature

# Push new branch to GitHub
git push -u origin feature/new-feature
```

#### Viewing Changes
```bash
# View file differences
git diff filename.py

# View staged changes
git diff --cached

# View changes between commits
git diff HEAD~1 HEAD

# View specific commit
git show commit-hash

# View file history
git log --follow filename.py
```

### GitHub Integration

#### Deploy Keys (Already Configured)
The server has SSH deploy keys configured for automatic deployment:
- **Public key**: Added to GitHub repository deploy keys
- **Private key**: Located at `/root/.ssh/repo_deploy`
- **SSH config**: Configured in `/root/.ssh/config`

#### GitHub Actions
The repository includes automated workflows:
- **Disk Space Monitor**: Runs every 6 hours
- **Location**: `.github/workflows/disk-space-monitor.yml`

#### Repository Secrets (For GitHub Actions)
Configure these in GitHub Settings → Secrets:
- `SERVER_HOST`: `api.moleculens.com`
- `SERVER_USER`: `root`
- `SERVER_SSH_KEY`: Private SSH key content

## 🚀 Deployment Operations

### Safe Deployment (Recommended)
```bash
# Use the safe deployment script
./scripts/maintenance/deploy-safe.sh
```

**What it does:**
1. Checks disk space before deployment
2. Runs cleanup if needed
3. Pulls latest code
4. Rebuilds containers safely
5. Performs health checks
6. Cleans up build artifacts

### Manual Deployment Steps
```bash
# 1. Kill processes on port 8000
fuser -k 8000/tcp || true

# 2. Navigate to project
cd /opt/hackathon-server/sci-vis-ai-server

# 3. Pull latest changes
git fetch origin main
git reset --hard origin/main

# 4. Navigate to Docker context
cd api

# 5. Stop containers
docker-compose down --remove-orphans

# 6. Clean up
docker container prune -f

# 7. Rebuild
docker-compose build

# 8. Start services
docker-compose up -d
```

### Rollback Deployment
```bash
# View recent commits
git log --oneline -10

# Rollback to specific commit
git reset --hard commit-hash

# Redeploy
./scripts/maintenance/deploy-safe.sh
```

## 🐳 Docker Operations

### Container Management
```bash
# List running containers
docker ps

# List all containers (including stopped)
docker ps -a

# View container logs
docker logs fastapi_app
docker logs --tail=100 fastapi_app
docker logs -f fastapi_app  # Follow logs

# Execute commands in container
docker exec -it fastapi_app bash
docker exec fastapi_app python -c "import sys; print(sys.version)"

# Restart container
docker restart fastapi_app

# Stop/start containers
docker-compose down
docker-compose up -d
```

### Image Management
```bash
# List images
docker images

# Remove unused images
docker image prune -f

# Remove specific image
docker rmi image-name

# Build image
docker-compose build

# Pull base images
docker pull continuumio/miniconda3
```

### Resource Monitoring
```bash
# View container resource usage
docker stats

# View container resource usage (one-time)
docker stats --no-stream

# Check disk usage by Docker
docker system df

# Detailed disk usage
docker system df -v
```

### Docker Cleanup
```bash
# Use automated cleanup script
./scripts/maintenance/docker-cleanup.sh

# Manual cleanup options
docker system prune -f                    # Remove unused data
docker system prune -a -f                 # Remove all unused data
docker system prune -a -f --volumes       # Remove everything unused

# Clean specific resources
docker container prune -f                 # Remove stopped containers
docker image prune -a -f                  # Remove unused images
docker volume prune -f                     # Remove unused volumes
docker network prune -f                   # Remove unused networks
```

## ⚙️ Service Management

### Application Service
```bash
# Check if application is running
curl http://localhost:8000/health

# Check application status
docker ps | grep fastapi_app

# View application logs
docker logs fastapi_app --tail=50

# Restart application
docker restart fastapi_app
```

### System Services
```bash
# Check Docker service
systemctl status docker

# Restart Docker service
systemctl restart docker

# Check SSH service
systemctl status ssh

# Check cron service
systemctl status cron

# View all services
systemctl list-units --type=service
```

### Maintenance Services
```bash
# Check maintenance timer
systemctl status moleculens-monitor.timer

# View maintenance logs
tail -f /opt/hackathon-server/sci-vis-ai-server/logs/maintenance/monitor.log

# Check cron jobs
crontab -l
```

## 📊 Monitoring & Debugging

### System Monitoring
```bash
# Check disk space
df -h

# Check memory usage
free -h

# Check CPU usage
top
htop  # If available

# Check system load
uptime

# Check running processes
ps aux | grep python
ps aux | grep docker
```

### Application Monitoring
```bash
# Test API endpoints
curl http://localhost:8000/health
curl http://localhost:8000/docs

# Check application performance
docker stats fastapi_app

# Monitor application logs in real-time
docker logs -f fastapi_app

# Check for errors in logs
docker logs fastapi_app 2>&1 | grep -i error
```

### Network Monitoring
```bash
# Check open ports
netstat -tlnp
ss -tlnp

# Check if port 8000 is in use
lsof -i :8000

# Test connectivity
ping google.com
curl -I http://localhost:8000
```

### Log Analysis
```bash
# Application logs
docker logs fastapi_app --since="1h"
docker logs fastapi_app --since="2023-08-09T10:00:00"

# System logs
journalctl -u docker --since="1 hour ago"
journalctl -f  # Follow system logs

# Maintenance logs
tail -f logs/maintenance/cleanup.log
tail -f logs/maintenance/monitor.log

# Search logs for errors
grep -i error logs/maintenance/*.log
journalctl -p err --since="1 day ago"
```

## 📁 File System Operations

### Navigation & File Management
```bash
# Navigate to project
cd /opt/hackathon-server/sci-vis-ai-server

# List files with details
ls -la

# Find files
find . -name "*.py" -type f
find . -name "*.log" -type f -mtime -1  # Modified in last day

# Check file sizes
du -sh *
du -sh logs/

# View file contents
cat filename.txt
less filename.txt  # For large files
head -20 filename.txt
tail -20 filename.txt
```

### File Editing
```bash
# Edit files with nano (beginner-friendly)
nano filename.py

# Edit files with vim (advanced)
vim filename.py

# Create/edit files with cat
cat > newfile.txt << EOF
Content here
EOF
```

### File Permissions
```bash
# Make file executable
chmod +x script.sh

# Set file permissions
chmod 644 file.txt      # Read/write for owner, read for others
chmod 755 script.sh     # Executable for owner, read/execute for others

# Change ownership
chown root:root file.txt

# View permissions
ls -la filename
```

### Backup Operations
```bash
# Backup important files
cp -r /opt/hackathon-server/sci-vis-ai-server /backup/sci-vis-ai-server-$(date +%Y%m%d)

# Create archive
tar -czf backup-$(date +%Y%m%d).tar.gz /opt/hackathon-server/sci-vis-ai-server

# Extract archive
tar -xzf backup.tar.gz
```

## 🔐 Security & Access

### SSH Key Management
```bash
# View SSH keys
ls -la ~/.ssh/

# Generate new SSH key
ssh-keygen -t ed25519 -C "your-email@example.com"

# View public key
cat ~/.ssh/id_ed25519.pub

# Test SSH connection
ssh -T git@github.com
```

### User Management
```bash
# Check current user
whoami
id

# View logged in users
who
w

# Check sudo access
sudo -l
```

### Firewall & Security
```bash
# Check firewall status
ufw status

# Check open ports
netstat -tlnp

# Check failed login attempts
grep "Failed password" /var/log/auth.log

# Check system security updates
apt list --upgradable
```

## 🚨 Emergency Procedures

### Service Recovery
```bash
# If application is down
docker restart fastapi_app

# If Docker is unresponsive
systemctl restart docker

# If system is unresponsive
reboot  # Use with caution
```

### Disk Space Emergency
```bash
# Check disk usage
df -h

# Emergency cleanup
./scripts/maintenance/docker-cleanup.sh --force

# Nuclear option (removes all unused Docker data)
docker system prune -a -f --volumes

# Clean system logs
journalctl --vacuum-time=7d

# Clean package cache
apt-get clean
```

### Git Recovery
```bash
# Recover from bad commit
git reflog
git reset --hard HEAD~1

# Recover deleted files
git checkout HEAD -- filename.py

# Force pull (loses local changes)
git fetch origin main
git reset --hard origin/main
```

### Application Recovery
```bash
# Reset to known good state
git reset --hard origin/main
./scripts/maintenance/deploy-safe.sh

# Check application health
curl http://localhost:8000/health

# View recent errors
docker logs fastapi_app --tail=100 | grep -i error
```

## 📚 Best Practices

### Development Workflow
1. **Always pull before making changes**:
   ```bash
   git pull origin main
   ```

2. **Use feature branches for major changes**:
   ```bash
   git checkout -b feature/new-feature
   ```

3. **Test locally before pushing**:
   ```bash
   docker-compose build
   docker-compose up -d
   ```

4. **Use descriptive commit messages**:
   ```bash
   git commit -m "feat: add user authentication endpoint"
   git commit -m "fix: resolve memory leak in image processing"
   git commit -m "docs: update API documentation"
   ```

### Deployment Best Practices
1. **Always use safe deployment**:
   ```bash
   ./scripts/maintenance/deploy-safe.sh
   ```

2. **Monitor deployments**:
   ```bash
   docker logs -f fastapi_app
   ```

3. **Verify health after deployment**:
   ```bash
   curl http://localhost:8000/health
   ```

4. **Keep backups of working versions**:
   ```bash
   git tag v1.0.0  # Tag stable versions
   ```

### Monitoring Best Practices
1. **Check disk space regularly**:
   ```bash
   df -h
   ```

2. **Monitor application logs**:
   ```bash
   tail -f logs/maintenance/monitor.log
   ```

3. **Review maintenance logs weekly**:
   ```bash
   tail -50 logs/maintenance/cleanup.log
   ```

### Security Best Practices
1. **Keep SSH keys secure**
2. **Use strong passwords**
3. **Regularly update system packages**:
   ```bash
   apt update && apt upgrade
   ```
4. **Monitor access logs**:
   ```bash
   tail -f /var/log/auth.log
   ```

## 🔍 Troubleshooting Common Issues

### Application Won't Start
```bash
# Check if port is in use
lsof -i :8000

# Check Docker status
docker ps
systemctl status docker

# Check application logs
docker logs fastapi_app

# Restart everything
docker-compose down
docker-compose up -d
```

### Git Issues
```bash
# Permission denied (publickey)
ssh -T git@github.com
cat ~/.ssh/config

# Merge conflicts
git status
git diff
# Edit conflicted files, then:
git add .
git commit -m "resolve merge conflicts"

# Detached HEAD state
git checkout main
git pull origin main
```

### Docker Issues
```bash
# Container won't start
docker logs container-name
docker inspect container-name

# Out of disk space
./scripts/maintenance/docker-cleanup.sh --force

# Permission issues
docker exec -it container-name ls -la /app
```

### Performance Issues
```bash
# Check resource usage
docker stats
top
df -h

# Check for memory leaks
docker logs fastapi_app | grep -i memory

# Restart services
docker restart fastapi_app
```

## 📞 Getting Help

### Log Locations
- **Application logs**: `docker logs fastapi_app`
- **Maintenance logs**: `logs/maintenance/*.log`
- **System logs**: `/var/log/syslog`, `journalctl`
- **Docker logs**: `journalctl -u docker`

### Useful Commands for Support
```bash
# System information
uname -a
df -h
free -h
docker --version
git --version

# Application status
docker ps
curl http://localhost:8000/health
git status
```

### Emergency Contacts
- **Repository**: [GitHub Issues](https://github.com/jakekinchen/moleculens-server/issues)
- **Documentation**: `SERVER_MAINTENANCE.md`
- **Server logs**: Check logs first before asking for help

---

## 🎯 Quick Reference

### Daily Commands
```bash
# Check system health
df -h && docker ps && curl -s http://localhost:8000/health

# Deploy safely
./scripts/maintenance/deploy-safe.sh

# View recent logs
docker logs --tail=20 fastapi_app
```

### Weekly Commands
```bash
# Review maintenance logs
tail -50 logs/maintenance/cleanup.log

# Check for updates
git log --oneline -10

# System cleanup (if needed)
./scripts/maintenance/docker-cleanup.sh --dry-run
```

### Emergency Commands
```bash
# Emergency cleanup
./scripts/maintenance/docker-cleanup.sh --force

# Reset to working state
git reset --hard origin/main && ./scripts/maintenance/deploy-safe.sh

# Restart everything
docker-compose down && docker-compose up -d
```

This guide covers all essential operations for managing the Moleculens API server. Keep this document handy for reference and update it as new procedures are added.