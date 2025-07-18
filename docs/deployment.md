# Deployment Guide

This guide explains how to deploy the Moleculens server to production.

## Prerequisites

- SSH access to the production server
- Docker and Docker Compose installed on the production server
- Domain name configured with DNS records
- SSL certificates (recommended)

## Deployment Steps

### 1. Initial Server Setup

```bash
# SSH into the server
ssh root@moleculens.com

# Create project directory
mkdir -p /opt/moleculens-server
cd /opt/moleculens-server

# Clone the repository
git clone https://github.com/yourusername/moleculens-server.git .
```

### 2. Environment Configuration

1. Create `.env` file:
```bash
cp .env.example .env
nano .env
```

Required environment variables:
```env
OPENAI_API_KEY=your-api-key
ENVIRONMENT=production
REDIS_HOST=redis
REDIS_PORT=6379
```

### 3. Docker Setup

1. Build and start the containers:
```bash
docker-compose up --build -d
```

2. Verify the containers are running:
```bash
docker-compose ps
```

### 4. Nginx Configuration (Recommended)

1. Install Nginx:
```bash
apt update
apt install nginx
```

2. Create Nginx configuration:
```nginx
server {
    listen 80;
    server_name moleculens.com www.moleculens.com;

    location / {
        proxy_pass http://localhost:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
    }
}
```

3. Enable the site and restart Nginx:
```bash
ln -s /etc/nginx/sites-available/moleculens /etc/nginx/sites-enabled/
nginx -t
systemctl restart nginx
```

### 5. SSL Configuration (Recommended)

1. Install Certbot:
```bash
apt install certbot python3-certbot-nginx
```

2. Obtain SSL certificate:
```bash
certbot --nginx -d moleculens.com -d www.moleculens.com
```

### 6. Maintenance and Updates

To update the server with latest changes:

```bash
# Pull latest changes
git pull origin main

# Rebuild and restart containers
docker-compose down
docker-compose up --build -d
```

### 7. Monitoring

1. Check container logs:
```bash
docker-compose logs -f
```

2. Monitor system resources:
```bash
docker stats
```

## Troubleshooting

### Common Issues

1. **Container fails to start**
   - Check logs: `docker-compose logs api`
   - Verify environment variables
   - Check disk space: `df -h`

2. **Rate limiting issues**
   - Check Redis connection
   - Verify rate limit settings in `api/utils/rate_limit.py`

3. **CORS errors**
   - Verify domain configuration in `api/main.py`
   - Check Nginx headers configuration

## Backup

1. Database backups (if applicable):
```bash
# Example backup command
docker-compose exec db pg_dump -U postgres > backup.sql
```

2. Environment configuration:
```bash
cp .env .env.backup
```

## Security Considerations

1. Keep the `.env` file secure and never commit it to version control
2. Regularly update dependencies
3. Monitor server logs for suspicious activity
4. Use strong passwords and API keys
5. Keep the system and Docker updated

## Support

For issues or questions, please:
1. Check the troubleshooting guide above
2. Review the logs
3. Open an issue on GitHub
4. Contact the development team
