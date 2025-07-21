"""PyMOL scene template helpers.

These functions return *deterministic* PyMOL command sequences for frequently
used visualisation tasks.  They **only** build strings – the caller is
responsible for executing them (e.g. via ``pymol.cmd``) inside a thread-safe
context.

Available helpers
-----------------
* ``overview_scene(structure_id)``
* ``binding_site_scene(structure_id, selection)``
* ``mutation_scene(structure_id, mutation_selection, original_residue=None)``
"""
from __future__ import annotations

from typing import List


# --------------------------------------------------------------------------- #
# Core scenes                                                                 #
# --------------------------------------------------------------------------- #

def overview_scene(structure_id: str) -> List[str]:
    """Generic overview – coloured cartoon plus transparent surface."""
    return [
        f"fetch {structure_id}, async=0",  # synchronous so caller may continue
        "hide everything",
        "show cartoon",
        # rainbow by chain for quick distinction
        "color chainbow, all",
        # semi-transparent molecular surface
        "set transparency, 0.4",
        "show surface",
        # white background & centred view
        "bg_color white",
        "orient",
    ]


def binding_site_scene(structure_id: str, selection: str) -> List[str]:
    """Highlight a binding-site selection.

    *structure_id*   PDB ID or object already present in PyMOL.
    *selection*      Valid PyMOL atom-selection string (e.g. ``\"resi 45+78 and chain A\"``).
    """
    return [
        f"fetch {structure_id}, async=0",
        "hide everything",
        "show cartoon",
        "color lightblue, all",
        # named selection for clarity
        "select binding_site, (" + selection + ")",
        "show sticks, binding_site",
        "color orange, binding_site",
        # focus the camera
        "set transparency, 0.5",
        "show surface, binding_site around 4",
        "zoom binding_site, 10",
        "bg_color white",
    ]


def mutation_scene(
    structure_id: str,
    mutation_selection: str,
    original_residue: str | None = None,
) -> List[str]:
    """Visualise a specific mutation site.

    Parameters
    ----------
    structure_id : str
        PDB ID or local object name to load.
    mutation_selection : str
        PyMOL selection string identifying the mutated residue(s),
        e.g. ``\"resi 123 and chain A\"``.
    original_residue : str, optional
        The one-letter or three-letter code of the wild-type residue,
        used for colouring if provided.

    Returns
    -------
    List[str]
        Ordered PyMOL commands.
    """
    colour_mut = "magenta"
    colour_orig = "yellow"

    cmds: List[str] = [
        f"fetch {structure_id}, async=0",
        "hide everything",
        "show cartoon",
        "color grey80, all",
        # highlight mutation
        f"select mutation_site, ({mutation_selection})",
        "show sticks, mutation_site",
        f"color {colour_mut}, mutation_site",
    ]

    if original_residue:
        # Attempt to colour wild-type residue differently if user supplied it
        cmds += [
            f"select wt_site, ({mutation_selection} and resn {original_residue})",
            f"color {colour_orig}, wt_site",
        ]

    cmds += [
        "set transparency, 0.35",
        "show surface, mutation_site around 4",
        "zoom mutation_site, 12",
        "bg_color white",
    ]
    return cmds


def mutation_focus_scene(structure_id: str, mutation_selection: str) -> List[str]:
    """Close-up view of a mutation site without surface.

    Parameters
    ----------
    structure_id : str
        PDB ID or local object name to load.
    mutation_selection : str
        PyMOL selection string targeting the mutated residue(s).

    Returns
    -------
    List[str]
        Ordered PyMOL commands.
    """

    return [
        f"fetch {structure_id}, async=0",
        "hide everything",
        "show cartoon",
        "color grey80, all",
        f"select mutation_site, ({mutation_selection})",
        "show sticks, mutation_site",
        "color magenta, mutation_site",
        "zoom mutation_site, 8",
        "bg_color white",
    ]
