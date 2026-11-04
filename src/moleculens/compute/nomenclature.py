"""Guarded Openclatura adapter for naming explanations."""

from __future__ import annotations

import os
import tempfile
from dataclasses import dataclass
from functools import wraps
from importlib.metadata import version
from pathlib import Path
from threading import Lock
from typing import Any, Literal

_UNSUPPORTED_METAL_ATOMIC_NUMBERS = frozenset(
    [
        3,
        4,
        11,
        12,
        13,
        19,
        20,
        31,
        37,
        38,
        49,
        50,
        55,
        56,
        81,
        82,
        83,
        84,
        87,
        88,
        113,
        114,
        115,
        116,
        *range(21, 31),
        *range(39, 49),
        *range(57, 81),
        *range(89, 113),
    ]
)
_MAX_ATOMS = 256
_PY2OPSIN_PATCH_LOCK = Lock()
OpsinStatus = Literal[
    "matched",
    "mismatched",
    "name_unparseable",
    "name_empty",
    "skipped_no_opsin",
    "skipped_no_java",
    "error",
    "not_requested",
]


def _safe_py2opsin(delegate: Any) -> Any:
    """Wrap py2opsin's shared cwd temp file with one unique file per call."""

    @wraps(delegate)
    def wrapped(*args: Any, **kwargs: Any) -> Any:
        if len(args) >= 8 or "tmp_fpath" in kwargs:
            return delegate(*args, **kwargs)
        descriptor, temp_path = tempfile.mkstemp(prefix="moleculens_opsin_", suffix=".txt")
        os.close(descriptor)
        try:
            return delegate(*args, tmp_fpath=temp_path, **kwargs)
        finally:
            Path(temp_path).unlink(missing_ok=True)

    setattr(wrapped, "__moleculens_safe_tmp__", True)  # noqa: B010
    return wrapped


def _ensure_safe_opsin_tempfiles() -> None:
    """Patch py2opsin once so concurrent API calls never share a cwd-relative file."""

    try:
        import py2opsin
    except ImportError:
        return
    with _PY2OPSIN_PATCH_LOCK:
        if getattr(py2opsin.py2opsin, "__moleculens_safe_tmp__", False):
            return
        py2opsin.py2opsin = _safe_py2opsin(py2opsin.py2opsin)


class NomenclatureError(Exception):
    """Base error for nomenclature requests."""


class UnsupportedNomenclatureError(NomenclatureError):
    """Raised when a molecule is outside the guarded initial domain."""

    def __init__(self, message: str, code: str) -> None:
        super().__init__(message)
        self.code = code


@dataclass(frozen=True)
class AtomMapping:
    """Mapping from Openclatura/RDKit input indices to viewer molblock indices."""

    status: Literal["matched", "unmatched", "not_requested"]
    engine_to_viewer: dict[int, int]


@dataclass(frozen=True)
class NomenclatureTrace:
    """One decision with separate engine and viewer atom domains."""

    phase: str
    decision: str
    reason: str
    atoms: list[int]
    engine_atoms: list[int]
    data: dict[str, Any]


@dataclass(frozen=True)
class NomenclatureExplanationComponent:
    """One human-readable explanation line."""

    phase: str
    text: str


@dataclass(frozen=True)
class NomenclatureExplanation:
    """Domain result returned by the Openclatura adapter."""

    display_name: str
    verified: bool
    opsin_status: OpsinStatus
    rules_hit: list[str]
    rule_hints: list[str]
    trace: list[NomenclatureTrace]
    components: list[NomenclatureExplanationComponent]
    atom_locants: dict[str, int]
    mapping_status: Literal["matched", "unmatched", "not_requested"]
    engine_version: str
    rdkit_version: str
    source: Literal["openclatura"] = "openclatura"
    experimental: bool = True


def _get_dependencies() -> tuple[Any, Any, Any, Any, Any]:
    _ensure_safe_opsin_tempfiles()
    from openclatura import describe, name
    from rdkit import Chem, rdBase

    return Chem, rdBase, name, describe, version


def _parse_and_validate_smiles(chem: Any, smiles: str) -> Any:
    mol = chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        raise NomenclatureError("Invalid SMILES")
    if mol.GetNumAtoms() > _MAX_ATOMS:
        raise UnsupportedNomenclatureError(
            f"Molecules with more than {_MAX_ATOMS} explicit atoms are not supported",
            "too_many_atoms",
        )
    fragments = chem.GetMolFrags(mol)
    if len(fragments) != 1:
        raise UnsupportedNomenclatureError(
            "Disconnected structures and mixtures are not supported yet",
            "disconnected_structure",
        )
    if any(atom.GetAtomicNum() in _UNSUPPORTED_METAL_ATOMIC_NUMBERS for atom in mol.GetAtoms()):
        raise UnsupportedNomenclatureError(
            "Metal-containing and coordination structures are not supported yet",
            "metal_structure",
        )
    return mol


def _map_to_viewer_atoms(chem: Any, engine_mol: Any, molblock: str | None) -> AtomMapping:
    if not molblock:
        return AtomMapping(status="not_requested", engine_to_viewer={})

    viewer_mol = chem.MolFromMolBlock(molblock, sanitize=True, removeHs=False)
    if viewer_mol is None:
        return AtomMapping(status="unmatched", engine_to_viewer={})

    engine_graph = chem.RemoveHs(engine_mol)
    viewer_graph = chem.RemoveHs(viewer_mol)
    if (
        engine_graph.GetNumAtoms() != viewer_graph.GetNumAtoms()
        or engine_graph.GetNumBonds() != viewer_graph.GetNumBonds()
        or chem.MolToSmiles(engine_graph, canonical=True, isomericSmiles=True)
        != chem.MolToSmiles(viewer_graph, canonical=True, isomericSmiles=True)
    ):
        return AtomMapping(status="unmatched", engine_to_viewer={})

    matches = viewer_mol.GetSubstructMatches(engine_mol, useChirality=True, uniquify=True)
    if not matches:
        matches = viewer_mol.GetSubstructMatches(engine_mol, useChirality=False, uniquify=True)
    if not matches:
        return AtomMapping(status="unmatched", engine_to_viewer={})

    match = min(matches)
    mapping = dict(enumerate(match))
    if len(mapping) != engine_mol.GetNumAtoms():
        return AtomMapping(status="unmatched", engine_to_viewer={})
    return AtomMapping(status="matched", engine_to_viewer=mapping)


def _viewer_atoms(atoms: Any, mapping: AtomMapping) -> list[int]:
    if mapping.status != "matched" or not isinstance(atoms, (list, tuple, set)):
        return []
    viewer_atoms = {
        mapping.engine_to_viewer[atom]
        for atom in atoms
        if isinstance(atom, int) and atom in mapping.engine_to_viewer
    }
    return sorted(viewer_atoms)


def _json_safe_trace_data(data: Any) -> dict[str, Any]:
    if not isinstance(data, dict):
        return {}
    safe: dict[str, Any] = {}
    for key, value in data.items():
        if isinstance(value, (str, int, float, bool)) or value is None:
            safe[str(key)] = value
        elif isinstance(value, (list, tuple)):
            safe[str(key)] = list(value)
        elif isinstance(value, dict):
            safe[str(key)] = {str(k): v for k, v in value.items()}
    return safe


def _build_trace(result: Any, mapping: AtomMapping) -> list[NomenclatureTrace]:
    trace: list[NomenclatureTrace] = []
    for step in result.decisions:
        engine_atoms = sorted(int(atom) for atom in step.atoms)
        trace.append(
            NomenclatureTrace(
                phase=step.phase.value,
                decision=step.decision,
                reason=step.reason,
                atoms=_viewer_atoms(engine_atoms, mapping),
                engine_atoms=engine_atoms,
                data=_json_safe_trace_data(step.data),
            )
        )
    return trace


def _walk_parent_nodes(nodes: Any) -> list[dict[str, Any]]:
    parents: list[dict[str, Any]] = []
    if not isinstance(nodes, (list, tuple)):
        return parents
    for node in nodes:
        if not isinstance(node, dict):
            continue
        parent = node.get("parent")
        if isinstance(parent, dict):
            parents.append(parent)
        for key in ("substituents", "replacement_prefixes", "instances"):
            parents.extend(_walk_parent_nodes(node.get(key)))
    return parents


def _extract_atom_locants(result: Any, mapping: AtomMapping) -> dict[str, int]:
    if mapping.status != "matched":
        return {}
    atom_locants: dict[str, int] = {}
    for parent in _walk_parent_nodes(result.substituent_tree):
        raw = parent.get("atom_ids_by_locant")
        if not isinstance(raw, dict):
            continue
        for locant, atom_id in raw.items():
            if not isinstance(atom_id, int) or atom_id not in mapping.engine_to_viewer:
                continue
            key = str(locant)
            if key not in atom_locants:
                atom_locants[key] = mapping.engine_to_viewer[atom_id]
    return atom_locants


def _format_locants(value: Any) -> str:
    if not isinstance(value, (list, tuple)):
        return ""
    return ", ".join(str(locant) for locant in value if str(locant))


def _build_components(result: Any, description: Any) -> list[NomenclatureExplanationComponent]:
    """Build educational prose without presenting engine atom ids as viewer ids."""

    components = [
        NomenclatureExplanationComponent(phase=component.phase, text=component.text)
        for component in description.components
        if component.phase not in {"parent_selection", "numbering"}
    ]
    bindings: list[dict[str, Any]] = []
    for step in result.decisions:
        raw_bindings = step.data.get("name_atom_bindings")
        if isinstance(raw_bindings, list):
            bindings.extend(binding for binding in raw_bindings if isinstance(binding, dict))

    for binding in bindings:
        stage = str(binding.get("stage") or "name part")
        role = str(binding.get("role") or "group")
        term = str(binding.get("term") or role)
        locants = _format_locants(binding.get("locants"))
        location = f" at {locants}" if locants else ""
        if stage == "parent":
            text = f"The {term} supplies the parent portion of the name."
        elif stage == "prefix" and role == "substituent":
            text = f"The {term} substituent{location} supplies the prefix."
        elif stage == "suffix":
            text = f"The principal {term} group{location} supplies the suffix."
        else:
            text = f"The {term} group{location} supplies the {stage} portion of the name."
        components.append(NomenclatureExplanationComponent(phase=stage, text=text))

    components.append(
        NomenclatureExplanationComponent(
            phase="assembly",
            text=f"The selected pieces assemble to form {result.name}.",
        )
    )
    return components


def explain_nomenclature(
    smiles: str,
    *,
    molblock: str | None = None,
    verify_opsin: bool = True,
) -> NomenclatureExplanation:
    """Generate a deterministic name explanation without leaking engine atom IDs as viewer IDs."""

    chem, rd_base, name_molecule, describe_molecule, package_version = _get_dependencies()
    _parse_and_validate_smiles(chem, smiles)

    result = name_molecule(smiles, include_trace=True, verify_opsin=verify_opsin)
    if not result.ok:
        raise NomenclatureError(result.error or "Openclatura did not produce a name")

    # Openclatura's atom ids belong to the SMILES stored on its result. Parse that
    # exact atom domain rather than assuming the caller's order survived processing.
    engine_mol = _parse_and_validate_smiles(chem, result.smiles)
    mapping = _map_to_viewer_atoms(chem, engine_mol, molblock)
    description = describe_molecule(result.smiles)
    opsin_status: OpsinStatus = (
        result.opsin_check.status if result.opsin_check is not None else "not_requested"
    )

    return NomenclatureExplanation(
        display_name=result.name,
        verified=result.verified,
        opsin_status=opsin_status,
        rules_hit=list(result.rules_hit),
        rule_hints=list(result.rule_hints),
        trace=_build_trace(result, mapping),
        components=_build_components(result, description),
        atom_locants=_extract_atom_locants(result, mapping),
        mapping_status=mapping.status,
        engine_version=package_version("openclatura"),
        rdkit_version=rd_base.rdkitVersion,
    )
