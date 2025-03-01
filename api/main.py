from fastapi import Depends, FastAPI
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
import os

# Try to load environment variables from .env file
try:
    from dotenv import load_dotenv
    load_dotenv()
    print("Environment variables loaded from .env file")
except ImportError:
    print("python-dotenv not installed. Environment variables must be set manually.")

import routers

app = FastAPI(
    title="AI Backend",
    version="0.0.1",
)
app.mount("/static", StaticFiles(directory="static"), name="static")

origins = [
    "*",
    "http://localhost:3000",
    "https://localhost:3000",
    "http://127.0.0.1:3000",
    "https://127.0.0.1:3000",
    "https://meshmo.com"
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["GET", "POST", "PUT", "DELETE", "OPTIONS"],
    allow_headers=["*"],
    expose_headers=["*"]
)



# user management related endpoints
app.include_router(routers.prompt.router)
app.include_router(routers.geometry.router)
