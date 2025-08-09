#!/usr/bin/env python3
"""
GitHub Webhook Receiver for Moleculens API
Listens for push events to main branch and triggers automatic deployment
"""

import os
import sys
import json
import hmac
import hashlib
import subprocess
import logging
from datetime import datetime
from http.server import HTTPServer, BaseHTTPRequestHandler
from urllib.parse import urlparse, parse_qs

# Configuration
WEBHOOK_SECRET = os.environ.get('WEBHOOK_SECRET', 'moleculens-webhook-secret-2024')
WEBHOOK_PORT = int(os.environ.get('WEBHOOK_PORT', '9000'))
PROJECT_DIR = '/opt/hackathon-server/sci-vis-ai-server'
DEPLOY_SCRIPT = f'{PROJECT_DIR}/scripts/maintenance/auto-deploy.sh'
LOG_FILE = f'{PROJECT_DIR}/logs/webhook.log'

# Setup logging
os.makedirs(os.path.dirname(LOG_FILE), exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

class WebhookHandler(BaseHTTPRequestHandler):
    def log_message(self, format, *args):
        """Override to use our logger"""
        logger.info(f"{self.address_string()} - {format % args}")

    def do_GET(self):
        """Handle GET requests (health check)"""
        if self.path == '/health':
            self.send_response(200)
            self.send_header('Content-type', 'application/json')
            self.end_headers()
            response = {
                'status': 'healthy',
                'service': 'moleculens-webhook',
                'timestamp': datetime.now().isoformat()
            }
            self.wfile.write(json.dumps(response).encode())
        else:
            self.send_response(404)
            self.end_headers()

    def do_POST(self):
        """Handle POST requests (GitHub webhooks)"""
        try:
            # Only accept webhooks on /webhook path
            if self.path != '/webhook':
                logger.warning(f"Invalid webhook path: {self.path}")
                self.send_response(404)
                self.end_headers()
                return

            # Get content length
            content_length = int(self.headers.get('Content-Length', 0))
            if content_length == 0:
                logger.warning("Empty webhook payload")
                self.send_response(400)
                self.end_headers()
                return

            # Read the payload
            payload = self.rfile.read(content_length)
            
            # Verify GitHub signature
            if not self.verify_signature(payload):
                logger.error("Invalid webhook signature")
                self.send_response(403)
                self.end_headers()
                return

            # Parse JSON payload
            try:
                data = json.loads(payload.decode('utf-8'))
            except json.JSONDecodeError as e:
                logger.error(f"Invalid JSON payload: {e}")
                self.send_response(400)
                self.end_headers()
                return

            # Process the webhook
            if self.process_webhook(data):
                self.send_response(200)
                self.send_header('Content-type', 'application/json')
                self.end_headers()
                response = {'status': 'success', 'message': 'Deployment triggered'}
                self.wfile.write(json.dumps(response).encode())
            else:
                self.send_response(200)
                self.send_header('Content-type', 'application/json')
                self.end_headers()
                response = {'status': 'ignored', 'message': 'No deployment needed'}
                self.wfile.write(json.dumps(response).encode())

        except Exception as e:
            logger.error(f"Error processing webhook: {e}")
            self.send_response(500)
            self.end_headers()

    def verify_signature(self, payload):
        """Verify GitHub webhook signature"""
        signature_header = self.headers.get('X-Hub-Signature-256')
        if not signature_header:
            logger.warning("Missing signature header")
            return False

        try:
            signature = signature_header.split('=')[1]
            expected_signature = hmac.new(
                WEBHOOK_SECRET.encode(),
                payload,
                hashlib.sha256
            ).hexdigest()
            
            return hmac.compare_digest(signature, expected_signature)
        except Exception as e:
            logger.error(f"Signature verification error: {e}")
            return False

    def process_webhook(self, data):
        """Process GitHub webhook data"""
        event_type = self.headers.get('X-GitHub-Event')
        
        logger.info(f"Received {event_type} event from GitHub")
        
        # Only process push events
        if event_type != 'push':
            logger.info(f"Ignoring {event_type} event")
            return False

        # Check if push is to main branch
        ref = data.get('ref', '')
        if ref != 'refs/heads/main':
            logger.info(f"Ignoring push to {ref} (not main branch)")
            return False

        # Get commit information
        commits = data.get('commits', [])
        if not commits:
            logger.info("No commits in push event")
            return False

        # Log commit details
        latest_commit = commits[-1]
        commit_hash = latest_commit.get('id', 'unknown')[:8]
        commit_message = latest_commit.get('message', 'No message')
        author = latest_commit.get('author', {}).get('name', 'Unknown')
        
        logger.info(f"Processing push to main: {commit_hash} by {author}")
        logger.info(f"Commit message: {commit_message}")

        # Trigger deployment
        return self.trigger_deployment(commit_hash, commit_message, author)

    def trigger_deployment(self, commit_hash, commit_message, author):
        """Trigger the deployment script"""
        try:
            logger.info("Starting automatic deployment...")
            
            # Check if deploy script exists
            if not os.path.exists(DEPLOY_SCRIPT):
                logger.error(f"Deploy script not found: {DEPLOY_SCRIPT}")
                return False

            # Set environment variables for the deployment
            env = os.environ.copy()
            env.update({
                'WEBHOOK_COMMIT_HASH': commit_hash,
                'WEBHOOK_COMMIT_MESSAGE': commit_message,
                'WEBHOOK_AUTHOR': author,
                'WEBHOOK_TIMESTAMP': datetime.now().isoformat()
            })

            # Run deployment script in background
            logger.info(f"Executing: {DEPLOY_SCRIPT}")
            
            # Use subprocess.Popen to run in background
            process = subprocess.Popen(
                [DEPLOY_SCRIPT],
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
                cwd=PROJECT_DIR
            )
            
            logger.info(f"Deployment started with PID: {process.pid}")
            
            # Don't wait for completion - let it run in background
            # The deployment script will handle its own logging
            
            return True
            
        except Exception as e:
            logger.error(f"Failed to trigger deployment: {e}")
            return False

def main():
    """Main function to start the webhook server"""
    logger.info(f"Starting Moleculens webhook receiver on port {WEBHOOK_PORT}")
    logger.info(f"Project directory: {PROJECT_DIR}")
    logger.info(f"Deploy script: {DEPLOY_SCRIPT}")
    
    # Check if deploy script exists
    if not os.path.exists(DEPLOY_SCRIPT):
        logger.error(f"Deploy script not found: {DEPLOY_SCRIPT}")
        logger.error("Please create the auto-deploy script first")
        sys.exit(1)
    
    try:
        server = HTTPServer(('0.0.0.0', WEBHOOK_PORT), WebhookHandler)
        logger.info(f"Webhook server ready at http://0.0.0.0:{WEBHOOK_PORT}")
        logger.info("Webhook endpoint: /webhook")
        logger.info("Health check: /health")
        server.serve_forever()
    except KeyboardInterrupt:
        logger.info("Shutting down webhook server...")
        server.shutdown()
    except Exception as e:
        logger.error(f"Server error: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()
