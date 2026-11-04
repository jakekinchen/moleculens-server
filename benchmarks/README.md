# Nomenclature benchmark

`nomenclature_corpus.json` is the rollout gate for Openclatura-backed naming. It combines
MolecuLens catalog and failure-prone structures with the round-trip cases distributed in the
locked MIT-licensed Openclatura package. The generated file records provenance per case and is
kept in git so benchmark inputs do not drift silently.

Refresh the corpus only when intentionally reviewing an engine upgrade:

```bash
uv run python scripts/build_nomenclature_corpus.py
```

Run a fast naming baseline without OPSIN:

```bash
uv run moleculens-nomenclature-benchmark benchmarks/nomenclature_corpus.json
```

Run the release gate in an environment with Java available:

```bash
uv run moleculens-nomenclature-benchmark benchmarks/nomenclature_corpus.json \
  --verify-opsin --enforce --output artifacts/nomenclature-benchmark.json
```

The enforced gate requires 300–500 cases, at least a 90% supported-name rate, 100% rejection of
the initial guarded domain, and at least a 95% OPSIN structural match rate. A passing OPSIN result
shows structural round-trip agreement; it does not prove that a generated string is the uniquely
preferred IUPAC name.
