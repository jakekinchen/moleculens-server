# AI Agent backend- fastapi, redis, and celery

## Overview
This repository demonstrates a Dockerized AI agent backend. 

## Prerequisites
- Docker and Docker Compose installed on your system.

## Setup Instructions

### 1. Create the `.env` File
Create a `.env` file in the root directory (same level as `docker-compose.yml`) and add the following variables:

```env
REDIS_HOST=redis
REDIS_PORT=6379
OPENAI_API_KEY=your-api-key
```

### 2. Build and Start the Containers
To build and start the containers, run:

```bash
docker-compose up --build
```
This command will build both the Redis and AI agent containers and start them.

### 3. Running the Agent
The AI agent will run and display its output in the terminal. Any code changes to `main.py` will reflect immediately without needing to rebuild.

### 4. Stopping the Containers
To stop and remove the containers, run:

```bash
docker-compose down
```

### 5. Rebuilding Without Cache
If you need to rebuild without using Docker's cache, use:

```bash
docker-compose build --no-cache
```

Now you're ready to build an AI agent


# To add the latest changes to the live server on meshmo.com do the following: 

## 1. ssh into the server 
```bash
    ssh root@meshmo.com
```

## 2. go to the project directory 
```bash 
    /opt/hackathon-server/sci-vis-ai-server/api
```

## 3. pull the latest changes 
```bash 
    git pull origin main
```

## 4. restart docker 
```bash
    docker-compose restart
```

