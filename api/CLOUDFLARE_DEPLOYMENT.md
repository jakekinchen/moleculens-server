# Deploying to Cloudflare Workers

This document explains how to deploy the Sci-Vis AI server to Cloudflare Workers using the Python runtime.

## Prerequisites

- Cloudflare account
- Wrangler CLI: `npm install -g wrangler`
- Node.js and npm

## Configuration Files

The project has been set up with the following Cloudflare Workers configuration files:

1. `wrangler.toml` - Main configuration file for Workers deployment
2. `workers_requirements.txt` - Documentation of package dependencies (note: only standard library is supported in production during beta)

## Environment Variables

The application requires several environment variables:

### API Keys (Sensitive)
- `OPENAI_API_KEY` - For OpenAI API access
- `ANTHROPIC_API_KEY` - For Anthropic/Claude API access
- `GROQ_API_KEY` - For Groq/Llama API access

### Celery Configuration
- `CELERY_BROKER_URL` - Redis URL for Celery broker (defaults to "redis://localhost:6379/0")
- `CELERY_RESULT_BACKEND` - Redis URL for Celery results (defaults to "redis://localhost:6379/0")

## Deployment Steps

1. **Set up the environment variables**

   For sensitive variables (API keys):
   ```
   npx wrangler secret put OPENAI_API_KEY
   npx wrangler secret put ANTHROPIC_API_KEY
   npx wrangler secret put GROQ_API_KEY
   ```

2. **Local testing**

   Run the application locally using Wrangler:
   ```
   npx wrangler dev
   ```

3. **Deploy to Cloudflare Workers**

   Deploy to Cloudflare Workers:
   ```
   npx wrangler deploy
   ```

## Limitations During Beta

During the Python Workers beta period:

1. Only standard library packages work in production
2. External packages (like FastAPI) only work in local development
3. The application has been adapted to handle both environments with conditional imports

## Celery Workers

The Celery workers are not deployed to Cloudflare Workers. They should be deployed separately as:

1. Standalone workers on a VPS/container platform
2. Adjusted to use a publicly accessible Redis instance or alternative message broker
3. Configured with the same API keys as the main application

## Monitoring and Logs

After deployment, you can view logs using:
```
npx wrangler tail
```

## Troubleshooting

If you encounter issues with the deployment:

1. Check Cloudflare Workers logs
2. Verify that all environment variables are set correctly
3. Ensure you're using the latest version of Wrangler
4. Check that you've added the `python_workers` compatibility flag in wrangler.toml