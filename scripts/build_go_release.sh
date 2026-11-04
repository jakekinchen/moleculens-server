#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
if [[ "$(go env GOVERSION)" != go1.27.3 ]]; then
  echo 'Release builder requires Go 1.27.3; update the pinned builder only with requalification.' >&2
  exit 1
fi
mkdir -p build
GOOS=linux GOARCH=amd64 CGO_ENABLED=0 go build -trimpath -o build/moleculens-linux ./cmd/moleculens
# Set DOCKER_CONTEXT for isolated local qualification, or leave Docker's selected context.
docker build --platform linux/amd64 -f docker/Dockerfile.engine -t moleculens-engine:qualified .
docker build --platform linux/amd64 -f docker/Dockerfile.go-api -t moleculens-go-api:candidate .
docker build --platform linux/amd64 -f docker/Dockerfile.go-worker -t moleculens-go-worker:candidate .
python3 - <<'PY'
import hashlib, json, subprocess
from pathlib import Path
images={}
for name in ('moleculens-engine:qualified','moleculens-go-api:candidate','moleculens-go-worker:candidate'):
    images[name]=subprocess.check_output(['docker','image','inspect',name,'--format','{{.Id}}'],text=True).strip()
receipt={'revision':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
         'dirty':bool(subprocess.check_output(['git','status','--porcelain'],text=True).strip()),
         'go':subprocess.check_output(['go','version'],text=True).strip(),
         'binarySha256':hashlib.sha256(Path('build/moleculens-linux').read_bytes()).hexdigest(),
         'images':images,'platform':'linux/amd64'}
Path('build/release.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
PY
