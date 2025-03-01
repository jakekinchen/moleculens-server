from fastapi import Depends, FastAPI
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware


import routers

app = FastAPI(
    title="AI Backend",
    version="0.0.1",
)
app.mount("/static", StaticFiles(directory="static"), name="static")

origins = [
    "*"
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)



# user management related endpoints
app.include_router(routers.prompt.router)
