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

import re
import tempfile
from typing import List

from api.agent_management.agents.pubchem_agent import PubChemAgent
from api.agent_management.llm_service import LLMService

# --------------------------------------------------------------------------- #
# Core scenes                                                                 #
# --------------------------------------------------------------------------- #

# Helper to determine if structure_id is a PDB ID
PDB_ID_PATTERN = re.compile(r"^[0-9a-zA-Z]{4}$")


def _get_structure_load_command(structure_id: str) -> str:
    """Return the appropriate PyMOL command to load a structure, supporting small molecules."""
    if PDB_ID_PATTERN.match(structure_id):
        return f"fetch {structure_id}, async=0"
    # Small molecule: fetch from PubChem, convert to PDB, save temp file, and load
    pubchem_agent = PubChemAgent(LLMService())
    try:
        pkg = pubchem_agent.get_molecule_package(structure_id)
        pdb_data = pkg.pdb_data
        if not pdb_data:
            raise ValueError(f"No PDB data for molecule: {structure_id}")
        with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp:
            tmp.write(pdb_data.encode("utf-8"))
            tmp_path = tmp.name
        return f"load {tmp_path}"
    except Exception as e:
        raise RuntimeError(f"Failed to fetch small molecule '{structure_id}': {e}")


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


def transparent_molecule_scene(structure_id: str, style: str = "cartoon") -> List[str]:
    """Generate a molecule with transparent background for overlays."""
    return [
        _get_structure_load_command(structure_id),
        "hide everything",
        f"show {style}",
        "color chainbow, all",
        "bg_color white",  # Will be made transparent by ray settings
        "orient",
        "set ray_opaque_background, 0",  # KEY: Enable transparency
        "set antialias, 1",
        "set ray_shadow, 1",
        "set depth_cue, 1",
    ]


def publication_quality_scene(
    structure_id: str, highlight_selection: str | None = None, transparent: bool = False
) -> List[str]:
    """High-quality rendering for publications with optional highlighting."""
    commands = [
        _get_structure_load_command(structure_id),
        "hide everything",
        "show cartoon",
        "color grey90, all",
        "set cartoon_highlight_color, blue",
        "set ray_trace_mode, 1",  # High quality mode
    ]
    if highlight_selection:
        commands.extend(
            [
                f"select highlight, ({highlight_selection})",
                "show sticks, highlight",
                "color hotpink, highlight",
                "set transparency, 0.3",
                "show surface, highlight around 4",
            ]
        )
    if transparent:
        commands.append("set ray_opaque_background, 0")
    return commands


def annotated_molecule_scene(
    structure_id: str, annotations: List[dict] | None = None
) -> List[str]:
    """Create molecule with labels, distances, and annotations."""
    commands = [
        _get_structure_load_command(structure_id),
        "hide everything",
        "show cartoon",
        "color chainbow, all",
        "orient",
        "set label_color, black",
        "set label_size, 14",
    ]
    if annotations:
        for annotation in annotations:
            if annotation["type"] == "distance":
                commands.append(
                    f"distance {annotation['name']}, {annotation['atom1']}, {annotation['atom2']}"
                )
            elif annotation["type"] == "label":
                commands.append(
                    f"label {annotation['selection']}, '{annotation['text']}'"
                )
            elif annotation["type"] == "angle":
                commands.append(
                    f"angle {annotation['name']}, {annotation['atom1']}, {annotation['atom2']}, {annotation['atom3']}"
                )
    return commands


def transparent_binding_site_scene(
    structure_id: str, selection: str, transparent_bg: bool = True
) -> List[str]:
    """Binding site visualization with transparent background for presentations."""
    commands = [
        _get_structure_load_command(structure_id),
        "hide everything",
        "show cartoon",
        "color lightblue, all",
        "select binding_site, (" + selection + ")",
        "show sticks, binding_site",
        "color orange, binding_site",
        "set transparency, 0.5",
        "show surface, binding_site around 4",
        "zoom binding_site, 10",
        "bg_color white",
        "set antialias, 1",
        "set ray_shadow, 1",
    ]
    if transparent_bg:
        commands.append("set ray_opaque_background, 0")
    return commands
