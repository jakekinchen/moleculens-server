# Instructions for fixing Docker issues:

## Fixing OpenAI Module Issue

The error `ModuleNotFoundError: No module named 'openai'` occurs because even though openai is in requirements.txt, there might be installation issues in the Docker container.

Try these solutions:

1. **Option 1: Update Dockerfile to install requirements during build**

Create a new `Dockerfile.fix` with:
```Dockerfile
FROM python:3.9

WORKDIR /app

# Install git and requirements during BUILD phase
COPY requirements.txt .
RUN apt-get update && apt-get install -y git && \
    pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# Set environment variables
ENV PYTHONPATH=/app

# Command to run the application with auto-reload
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000", "--reload"]
```

Then rebuild and restart:
```bash
cd /Users/jakekinchen/Documents/sci-vis-ai-server/api
docker-compose build --no-cache api
docker-compose up -d
```

2. **Option 2: Install openai package directly in the container**
```bash
docker exec -it fastapi_app pip install openai
docker-compose restart api
```

3. **Option 3: Specify openai version in requirements.txt**
If the package version is causing issues, try specifying a version:
```
openai==1.3.7
```
Then rebuild the containers.