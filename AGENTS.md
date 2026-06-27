# MolecuLens Server Agent Guide

This file is the root instruction layer for agents working in `moleculens-server`, the FastAPI service that computes and serves molecular orbital/density overlays and conformer artifacts for the MolecuLens clients.

## Project Shape

- API: FastAPI routes under `src/moleculens/api/**`
- Worker: job processing under `src/moleculens/worker/**`
- Database/job queue: SQLAlchemy/Postgres-backed async processing
- Cache/artifacts: deterministic geometry + parameter cache outputs
- Tests: `tests/unit/**` and `tests/integration/**`
- Tooling: `uv`, `ruff`, `mypy`, `pytest`, Docker Compose

## Non-Negotiables

- Do not change API payload shapes without checking client compatibility in the web/native repos.
- Cache keys must be deterministic and include all parameters that affect output.
- API code must not assume worker-local filesystem state unless the contract explicitly guarantees shared storage.
- Worker changes that affect queue, cache, or artifact behavior require integration proof.
- Do not commit secrets, production database URLs, metrics passwords, provider credentials, or local machine paths.
- Keep broad cleanup out of scoped bugfix branches.

## Golden Commands

```bash
uv sync --extra dev
uv run ruff format --check .
uv run ruff check .
uv run mypy src
uv run pytest tests/unit -v
```

Integration path:

```bash
docker compose -f docker-compose.yml -f docker-compose.dev.yml up -d
uv run pytest tests/integration -v
docker compose -f docker-compose.yml -f docker-compose.dev.yml down -v
```

## QA Matrix

| Change type | Required proof |
| --- | --- |
| Docs-only | Review only; state checks not run |
| Formatting/lint-only | `uv run ruff format --check .` + `uv run ruff check .` |
| Pure utility/model logic | ruff + mypy + targeted unit tests |
| API route behavior | ruff + mypy + unit tests + integration/API proof when behavior changes |
| Queue/worker/cache/artifacts | ruff + mypy + unit tests + Docker integration tests |
| Payload/client contract | server proof + client compatibility note for `jakekinchen/moleculens` |
| Docker/deploy config | Docker build/up proof or a clear reason it could not be run |

## Client Compatibility Checklist

Before changing orbital/conformer response contracts, verify or document:

- cache key fields,
- job status states,
- artifact URL format,
- mesh encoding/compression shape,
- metadata field names,
- error response shape,
- CORS/origin expectations,
- timeout/retry implications for web and native clients.

## PR Proof Bundle

Use this format in PR bodies:

```md
## Validation
- `uv run ruff format --check .`: pass/fail/not run
- `uv run ruff check .`: pass/fail/not run
- `uv run mypy src`: pass/fail/not run
- Unit tests: pass/fail/not run — <command>
- Integration tests: pass/fail/not run — <command>
- Client compatibility: not applicable / documented / tested

## Risk
- <what could regress>

## Follow-ups
- <deferred tasks>
```

Never claim a command passed unless it ran on the branch being submitted.

## Parallel Work

Use one branch/worktree per task. Treat these as mutex areas unless a sweep owner explicitly coordinates merge order:

- database migrations,
- Docker Compose/deploy files,
- cache key logic,
- public response models,
- worker queue semantics,
- dependency lock/config files.

## Completion Standard

A task is complete when:

1. the diff is scoped,
2. required proof is recorded,
3. client contract impact is documented,
4. risks/follow-ups are in the PR body,
5. any reusable rule is promoted back into this file or the client repo agent ops docs.
