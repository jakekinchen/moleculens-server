FROM continuumio/miniconda3

RUN conda install -y -c conda-forge python=3.9 pymol-open-source=3.1.0 assimp ffmpeg && conda clean -afy

WORKDIR /app
COPY api/requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

COPY . /app
RUN rm -f /app/__init__.py

CMD ["gunicorn", "api.main:app", "-k", "uvicorn.workers.UvicornWorker", "-b", "0.0.0.0:8000", "-w", "1"]
