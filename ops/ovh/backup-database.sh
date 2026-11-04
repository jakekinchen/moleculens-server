#!/usr/bin/env bash
set -euo pipefail
[[ "$(id -u)" == 0 ]] || { echo 'Run this backup with sudo.' >&2; exit 1; }
cd /opt/moleculens
umask 077
backup_dir=/var/backups/moleculens
install -d -m 0700 "$backup_dir"
exec 9>/run/moleculens-database-backup.lock
flock -n 9 || { echo 'A database backup is already running.' >&2; exit 1; }
available=$(df -B1 --output=avail "$backup_dir" | tail -n 1 | tr -d ' ')
(( available > 5368709120 )) || { echo 'Less than 5 GiB free: investigate storage before backing up.' >&2; exit 1; }
stamp=$(date -u +%Y%m%dT%H%M%SZ)
partial="$backup_dir/db-$stamp.dump.partial"
final="$backup_dir/db-$stamp.dump"
dc=(docker compose --env-file .env -f compose.vps1.yml)
"${dc[@]}" exec -T postgres pg_dump -U moleculens -d moleculens --format=custom < /dev/null > "$partial"
"${dc[@]}" exec -T postgres pg_restore --list < "$partial" > /dev/null
mv "$partial" "$final"
sha256sum "$final" > "$final.sha256"
# Rotate only this script's completed dumps, keeping the seven newest.
python3 - "$backup_dir" <<'PY'
import re, sys
from pathlib import Path
root=Path(sys.argv[1])
files=sorted(p for p in root.iterdir() if re.fullmatch(r'db-\d{8}T\d{6}Z\.dump',p.name))
for old in files[:-7]:
    old.unlink()
    old.with_name(old.name+'.sha256').unlink(missing_ok=True)
PY
printf '%s\n' "$final"
