"""Deterministic prompt translator for PyMOL.

This module converts very simple text prompts into pre-defined PyMOL command
sequences provided by :mod:`pymol_templates`.  It serves as a lightweight
alternative to the LLM-based translator used elsewhere.
"""
from __future__ import annotations

from typing import List

from .pymol_templates import binding_site_scene, overview_scene, mutation_scene


def translate(prompt: str) -> List[str]:
    """Translate a short text prompt into PyMOL commands."""
    words = prompt.split()
    if not words:
        return []

    if words[0] == "overview" and len(words) >= 2:
        struct_id = words[1]
        return overview_scene(struct_id)

    if words[0] == "binding" and len(words) >= 4 and words[1] == "site":
        struct_id = words[2]
        selection_str = " ".join(words[3:])
        return binding_site_scene(struct_id, selection_str)

    # mutation-focus prompt, e.g. "mutation 1abc resi 123 and chain A"
    if words[0] in {"mutation", "mutate"} and len(words) >= 4:
        struct_id = words[1]
        selection_str = " ".join(words[2:])
        return mutation_scene(struct_id, selection_str)

    return []
