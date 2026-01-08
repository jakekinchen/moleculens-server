"""Conformer generation for canonical 2D molblocks."""

import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from moleculens.core import get_logger

logger = get_logger(__name__)


@dataclass
class ConformerResult:
    """Result of conformer computation."""

    sdf3d: str
    meta: dict[str, Any] = field(default_factory=dict)
    artifact_files: list[str] = field(default_factory=list)
    compute_time_ms: float = 0.0


class ConformerComputationError(Exception):
    """Raised when conformer computation fails."""


def _get_rdkit() -> tuple[Any, Any]:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as e:
        raise ConformerComputationError("RDKit not available in this environment") from e
    return Chem, AllChem


def _build_etkdg_params(all_chem: Any, method: str, max_attempts: int) -> tuple[Any, str]:
    method_norm = method.lower()
    if method_norm in ("etkdg_v3", "etkdgv3", "etkdg3"):
        params = all_chem.ETKDGv3()
        method_norm = "etkdg_v3"
    elif method_norm in ("etkdg_v2", "etkdgv2", "etkdg2"):
        params = all_chem.ETKDGv2()
        method_norm = "etkdg_v2"
    else:
        params = all_chem.ETKDG()
        method_norm = "etkdg"
    return params, method_norm


def _has_metal(chem: Any, mol: Any) -> bool:
    try:
        periodic_table = chem.GetPeriodicTable()
    except Exception:
        return False

    for atom in mol.GetAtoms():
        try:
            if periodic_table.IsMetal(atom.GetAtomicNum()):
                return True
        except Exception:
            return False
    return False


def _ensure_planar_coords(chem: Any, mol: Any) -> None:
    if mol.GetNumConformers() == 0:
        conf = chem.Conformer(mol.GetNumAtoms())
        for idx in range(mol.GetNumAtoms()):
            conf.SetAtomPosition(idx, (0.0, 0.0, 0.0))
        mol.AddConformer(conf, assignId=True)
        return

    conf = mol.GetConformer()
    for idx in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(idx)
        conf.SetAtomPosition(idx, (pos.x, pos.y, 0.0))


def _optimize_molecule(all_chem: Any, mol: Any, opt: str, max_opt_iters: int) -> str:
    opt_norm = opt.lower()
    if opt_norm == "uff":
        if not all_chem.UFFHasAllMoleculeParams(mol):
            return "params_missing"
        result = all_chem.UFFOptimizeMolecule(mol, maxIters=max_opt_iters)
        return "optimized" if result == 0 else "not_converged"
    if opt_norm == "mmff":
        if not all_chem.MMFFHasAllMoleculeParams(mol):
            return "params_missing"
        result = all_chem.MMFFOptimizeMolecule(mol, maxIters=max_opt_iters)
        return "optimized" if result == 0 else "not_converged"
    return "skipped"


def run_conformer_computation(
    molblock2d: str,
    method: str = "etkdg_v3",
    opt: str = "uff",
    max_attempts: int = 10,
    max_opt_iters: int = 200,
    add_hs: bool = False,
    geom_version: str = "v1",
    output_dir: Path | None = None,
) -> ConformerResult:
    """Generate a 3D conformer from a 2D molblock.

    The output preserves atom order from the input molblock.
    """
    start_time = time.time()
    timings: dict[str, float] = {}

    Chem, AllChem = _get_rdkit()

    mol = Chem.MolFromMolBlock(molblock2d, sanitize=True, removeHs=False)
    if mol is None:
        raise ConformerComputationError("Failed to parse molblock2d")
    if mol.GetNumAtoms() == 0:
        raise ConformerComputationError("No atoms in molecule")

    has_metal = _has_metal(Chem, mol)

    opt_norm = opt.lower() if opt else "none"
    if opt_norm in ("off", "false", "no"):
        opt_norm = "none"

    params, method_norm = _build_etkdg_params(AllChem, method, max_attempts)

    logger.info(
        "Starting conformer computation",
        method=method_norm,
        opt=opt_norm,
        num_atoms=mol.GetNumAtoms(),
        add_hs=add_hs,
    )

    embed_status = "failed"
    opt_status = "skipped"
    quality = "planar_fallback"

    t0 = time.time()
    embed_result = 1
    for attempt in range(max_attempts):
        if attempt > 0:
            params.useRandomCoords = True
            mol.RemoveAllConformers()
        embed_result = AllChem.EmbedMolecule(mol, params)
        if embed_result == 0:
            break
    timings["embed"] = (time.time() - t0) * 1000

    if embed_result == 0:
        embed_status = "embedded"

        if add_hs:
            mol = Chem.AddHs(mol, addCoords=True)

        if opt_norm != "none":
            if has_metal:
                opt_status = "skipped_metal"
            else:
                t0 = time.time()
                opt_status = _optimize_molecule(AllChem, mol, opt_norm, max_opt_iters)
                timings["opt"] = (time.time() - t0) * 1000
        quality = "optimized" if opt_status == "optimized" else "embedded"
    else:
        _ensure_planar_coords(Chem, mol)
        if add_hs:
            mol = Chem.AddHs(mol, addCoords=True)
            _ensure_planar_coords(Chem, mol)

    sdf3d = Chem.MolToMolBlock(mol)

    total_time = (time.time() - start_time) * 1000
    timings["total"] = total_time

    meta = {
        "method": method_norm,
        "opt": opt_norm,
        "max_attempts": max_attempts,
        "max_opt_iters": max_opt_iters,
        "add_hs": add_hs,
        "geom_version": geom_version,
        "quality": quality,
        "has_metal": has_metal,
        "embed_status": embed_status,
        "opt_status": opt_status,
        "timings_ms": timings,
    }

    result = ConformerResult(sdf3d=sdf3d, meta=meta, compute_time_ms=total_time)

    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)
        structure_path = output_dir / "structure.sdf"
        structure_path.write_text(sdf3d)
        result.artifact_files.append("structure.sdf")

    logger.info(
        "Conformer computation complete",
        quality=quality,
        total_time_ms=total_time,
    )

    return result
