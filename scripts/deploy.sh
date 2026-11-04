#!/usr/bin/env bash
set -euo pipefail
echo 'The unfenced legacy deployment path is retired for this migration.' >&2
echo 'Build with scripts/build_go_release.sh and follow docs/go-release-runbook.md.' >&2
echo 'Verify the host identity and complete the recorded drain before any live cutover.' >&2
exit 2
