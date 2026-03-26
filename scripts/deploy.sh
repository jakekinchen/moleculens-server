#!/usr/bin/env bash

set -euo pipefail

compose_file="${DEPLOY_COMPOSE_FILE:-}"
if [[ -z "$compose_file" ]]; then
  if [[ -f docker-compose.single.yml ]]; then
    compose_file="docker-compose.single.yml"
  else
    compose_file="docker-compose.yml"
  fi
fi

echo "[deploy] using compose file ${compose_file}"

docker network create infra_app >/dev/null 2>&1 || true

docker compose -f "${compose_file}" --env-file .env build --progress plain
docker compose -f "${compose_file}" --env-file .env up -d --remove-orphans
docker compose -f "${compose_file}" --env-file .env ps

if docker ps --format '{{.Names}}' | grep -qx molecule-api; then
  echo "[deploy] waiting for api health"
  for i in 1 2 3 4 5 6 7 8 9 10 11 12; do
    if docker exec molecule-api curl --max-time 5 -sf http://localhost:8000/health >/dev/null 2>&1; then
      echo "[deploy] api healthy"
      docker image prune -f >/dev/null 2>&1 || true
      exit 0
    fi
    sleep 5
  done
  docker logs --tail=100 molecule-api || true
  exit 1
fi

echo "[deploy] molecule-api container not found"
docker ps -a
exit 1
