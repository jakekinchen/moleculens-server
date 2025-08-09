import hashlib
import json
import os
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional

from fastapi import APIRouter, HTTPException

router = APIRouter(prefix="/v1/figure", tags=["FigureSpec v1"])


# ----------------------------------------------------------------------------
# Deterministic canonicalization (MUST MATCH TS LIB)
# ----------------------------------------------------------------------------


def _round(x: Any):
    if isinstance(x, float):
        return round(x, 6)
    if isinstance(x, list):
        return [_round(i) for i in x]
    if isinstance(x, dict):
        return {k: _round(x[k]) for k in sorted(x)}
    return x


def canonical_json(obj: dict) -> str:
    return json.dumps(_round(obj), separators=(",", ":"), ensure_ascii=False)


def compute_spec_id(spec: dict) -> str:
    return hashlib.sha256(canonical_json(spec).encode("utf-8")).hexdigest()


# ----------------------------------------------------------------------------
# In-memory store and background processing
# ----------------------------------------------------------------------------


@dataclass
class FigureRecord:
    spec_id: str
    status: str = "queued"  # queued | processing | completed | failed
    urls: Optional[Dict[str, str]] = None
    error: Optional[str] = None


_records: Dict[str, FigureRecord] = {}
_records_lock = threading.Lock()


def _static_base_url() -> str:
    # Prefer explicit base URL; default to required production base
    return os.environ.get("PYMOL_SERVER_BASE_URL", "https://api.moleculens.com").rstrip("/")


def _asset_dir(spec_id: str) -> Path:
    # Assets live under api/static/figures/{spec_id}
    base = Path(__file__).resolve().parents[2] / "static" / "figures" / spec_id
    base.mkdir(parents=True, exist_ok=True)
    return base


def _asset_urls(spec_id: str) -> Dict[str, str]:
    base = _static_base_url()
    prefix = f"{base}/static/figures/{spec_id}"
    urls: Dict[str, str] = {}
    if (_asset_dir(spec_id) / "2d.svg").exists():
        urls["svg2d"] = f"{prefix}/2d.svg"
    if (_asset_dir(spec_id) / "2d.png").exists():
        urls["png2d"] = f"{prefix}/2d.png"
    if (_asset_dir(spec_id) / "3d.png").exists():
        urls["png3d"] = f"{prefix}/3d.png"
    if (_asset_dir(spec_id) / "scene.glb").exists():
        urls["glb"] = f"{prefix}/scene.glb"
    if (_asset_dir(spec_id) / "meta.json").exists():
        urls["meta"] = f"{prefix}/meta.json"
    return urls


def _write_once(path: Path, content: bytes) -> None:
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    # Atomic-ish write
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_bytes(content)
    tmp.replace(path)


def _generate_2d_molecule_svg(spec: dict, width: int, height: int, transparent: bool, spec_id: str) -> bytes:
    """Generate actual 2D molecular depiction using RDKit."""
    try:
        # Import RDKit drawing modules
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit.Chem.Draw import rdMolDraw2D

        # Extract input parameters
        input_data = spec.get("input", {})
        kind = input_data.get("kind", "smiles")
        value = input_data.get("value", "")

        # Parse molecule based on input type
        mol = None
        if kind == "smiles":
            mol = Chem.MolFromSmiles(value)
        elif kind == "name":
            # For molecule names, we'd need to resolve to SMILES first
            # For now, treat as SMILES if it looks like one, otherwise create placeholder
            mol = Chem.MolFromSmiles(value)
        elif kind == "pdb":
            # For PDB input, we'd need different parsing
            # For now, create placeholder
            pass

        if mol is None:
            # Fallback to placeholder if molecule parsing fails
            return _generate_placeholder_svg(width, height, transparent, spec_id, f"Invalid {kind}: {value}")

        # Generate 2D coordinates if not present
        if mol.GetNumConformers() == 0:
            AllChem.Compute2DCoords(mol)  # type: ignore[attr-defined]

        # Create SVG drawer
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)

        # Configure drawing options
        opts = drawer.drawOptions()
        if transparent:
            opts.clearBackground = False  # Don't draw background
        else:
            opts.backgroundColour = (1, 1, 1, 1)  # White background

        # Set font size based on image size
        opts.baseFontSize = max(12, min(24, width // 40))

        # Draw the molecule
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()

        # Get SVG content
        svg_content = drawer.GetDrawingText()
        return svg_content.encode("utf-8")

    except Exception as e:
        # Log error and return placeholder
        print(f"Error generating 2D molecule depiction: {e}")
        return _generate_placeholder_svg(width, height, transparent, spec_id, f"Rendering error: {str(e)}")


def _generate_placeholder_svg(width: int, height: int, transparent: bool, spec_id: str, error_msg: str = "") -> bytes:
    """Generate placeholder SVG when molecular rendering fails."""
    svg_bg = "none" if transparent else "white"
    display_msg = error_msg if error_msg else f"Moleculens Figure {spec_id[:8]}"

    svg = (
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>"
        f"<rect width='100%' height='100%' fill='{svg_bg}'/>"
        f"<text x='{width // 2}' y='{height // 2}' dominant-baseline='middle' text-anchor='middle' font-family='Helvetica' font-size='16'>"
        f"{display_msg}"
        f"</text></svg>"
    )
    return svg.encode("utf-8")


def _render_assets(spec: dict, spec_id: str) -> None:
    # Real molecular rendering pipeline replacing placeholders, honoring modes
    asset_dir = _asset_dir(spec_id)

    # Extract rendering parameters
    render_obj = spec.get("render", {}) or {}
    modes = [str(m).lower() for m in (render_obj.get("modes") or [])]
    want_2d = ("2d" in modes) or not modes  # default to 2D if modes not specified
    want_3d = "3d" in modes

    width = int(render_obj.get("width", 1024) or 1024)
    height = int(render_obj.get("height", 768) or 768)
    transparent = bool(render_obj.get("transparent", True))
    dpi = int(render_obj.get("dpi", 300) or 300)

    # Collect meta diagnostics
    errors: list[str] = []
    renderers: dict[str, str] = {}

    # Generate 2D molecular depiction only if requested
    if want_2d:
        svg_content = _generate_2d_molecule_svg(spec, width, height, transparent, spec_id)
        _write_once(asset_dir / "2d.svg", svg_content)

        # 2D PNG via cairosvg
        try:
            import cairosvg  # type: ignore

            png_bytes = cairosvg.svg2png(
                bytestring=svg_content,
                output_width=width,
                output_height=height,
                dpi=dpi,
            )
            _write_once(asset_dir / "2d.png", png_bytes)
            renderers["2d"] = "rdkit"
        except Exception as e:
            print(f"[{spec_id}] 2D PNG conversion failed: {e}")
            errors.append(f"2d_png_error: {e}")

    # 3D PNG generation only if requested; no silent 2D duplication
    if want_3d:
        three_d_ok = False
        try:
            # Build a PDB block from inputs
            input_obj = spec.get("input", {})
            kind = input_obj.get("kind")
            value = input_obj.get("value")
            # conformer_method is accepted at input-level but we always embed for 3D

            pdb_block: Optional[str] = None

            if kind == "smiles":
                try:
                    from rdkit import Chem  # type: ignore
                    from rdkit.Chem import AllChem  # type: ignore

                    mol = Chem.MolFromSmiles(value)
                    if mol is not None:
                        mol = Chem.AddHs(mol)
                        # For 3D rendering we must have 3D coordinates regardless of conformer_method
                        # Always attempt ETKDG embedding and MMFF optimization
                        try:
                            AllChem.EmbedMolecule(mol, AllChem.ETKDG())  # type: ignore[attr-defined]
                            try:
                                AllChem.MMFFOptimizeMolecule(mol)  # type: ignore[attr-defined]
                            except Exception:
                                pass
                        except Exception as e:
                            errors.append(f"rdkit_embed_failed: {e}")
                        if mol.GetNumConformers() > 0:
                            pdb_block = Chem.MolToPDBBlock(mol)
                        else:
                            errors.append("rdkit_conformer_failed")
                except Exception as e:
                    errors.append(f"rdkit_error: {e}")

            elif kind == "name" and isinstance(value, str) and value:
                # Resolve name via PubChem and fetch SDF, then convert to PDB
                try:
                    from api.pymol.services.pubchem import PubChemSearchService  # local import to avoid import cycle
                    from api.pymol.services.rdkit_utils import sdf_to_pdb_block  # local import

                    svc = PubChemSearchService()
                    compounds = svc.search_with_fallbacks(value)
                    if not compounds:
                        errors.append("pubchem_no_results")
                    else:
                        cid = getattr(compounds[0], "cid", None)
                        if not isinstance(cid, int):
                            errors.append("pubchem_invalid_cid")
                        else:
                            sdf_text = svc.fetch_sdf(cid, "3d") or svc.fetch_sdf(cid, "2d")
                            if not sdf_text:
                                errors.append("pubchem_no_sdf")
                            else:
                                pdb_from_sdf = sdf_to_pdb_block(sdf_text)
                                if pdb_from_sdf:
                                    pdb_block = pdb_from_sdf
                                else:
                                    errors.append("sdf_to_pdb_failed")
                except Exception as e:
                    errors.append(f"pubchem_error: {e}")

            elif kind == "pdb" and isinstance(value, str) and value:
                # Treat as PDB identifier and fetch from RCSB
                try:
                    import requests  # type: ignore

                    url = f"https://files.rcsb.org/download/{value.upper()}.pdb"
                    resp = requests.get(url, timeout=20)
                    if resp.ok and resp.text:
                        pdb_block = resp.text
                    else:
                        errors.append(f"rcsb_fetch_failed: {resp.status_code}")
                except Exception as e:
                    errors.append(f"rcsb_error: {e}")

            # If we have a PDB block, render with PyMOL
            if pdb_block:
                try:
                    import os as _os

                    _os.environ.setdefault("PYMOL_QUIET", "1")
                    _os.environ.setdefault("PYMOL_HEADLESS", "1")
                    import pymol  # type: ignore

                    if hasattr(pymol, "finish_launching"):
                        pymol.finish_launching(["pymol", "-cq"])  # quiet, headless
                    cmd = getattr(pymol, "cmd", None)
                    if cmd is None:
                        errors.append("pymol_cmd_missing")
                    else:
                        three_d = spec.get("3d", {}) or {}
                        representation = str(three_d.get("representation", "licorice")).lower()
                        bg = str(three_d.get("bg", "transparent")).lower()

                        cmd.reinitialize()
                        if hasattr(cmd, "read_pdbstr"):
                            cmd.read_pdbstr(pdb_block, "mol")
                        else:
                            tmp = asset_dir / "_tmp.pdb"
                            tmp.write_text(pdb_block)
                            cmd.load(str(tmp), "mol")
                            try:
                                tmp.unlink()
                            except Exception:
                                pass

                        cmd.hide("everything")
                        if "cartoon+licorice" in representation:
                            cmd.show("cartoon", "mol")
                            cmd.show("sticks", "mol")
                        elif "surface" in representation:
                            cmd.show("surface", "mol")
                        else:
                            cmd.show("sticks", "mol")
                        cmd.orient("mol")

                        # Background
                        if transparent:
                            cmd.set("ray_opaque_background", 0)
                        if bg == "black":
                            cmd.do("bg_color black")
                        elif bg == "white":
                            cmd.do("bg_color white")
                        else:
                            cmd.do("bg_color white")

                        try:
                            cmd.ray(width, height)
                        except Exception as e:
                            errors.append(f"pymol_ray_error: {e}")
                        out_path = asset_dir / "3d.png"
                        try:
                            cmd.png(str(out_path), dpi=dpi, ray=1)
                            three_d_ok = True
                            renderers["3d"] = "pymol"
                        except Exception as e:
                            errors.append(f"pymol_png_error: {e}")
                except Exception as e:
                    errors.append(f"pymol_error: {e}")
            else:
                errors.append("no_pdb_block")
        except Exception as e:
            errors.append(f"3d_pipeline_error: {e}")

        # Do not silently duplicate 2D into 3D. Only create 3d.png on success.
        if not three_d_ok:
            renderers["3d"] = "none"

    # Meta JSON (write-once)
    meta = {
        "spec_id": spec_id,
        "canonical": json.loads(canonical_json(spec)),
        "inputs": spec.get("input", {}),
        "render": spec.get("render", {}),
        "style_preset": spec.get("style_preset"),
        "annotations": spec.get("annotations", {}),
        "3d": spec.get("3d", {}),
        "renderers": renderers,
        "errors": errors,
    }
    _write_once(asset_dir / "meta.json", json.dumps(meta, separators=(",", ":")).encode("utf-8"))


def _validate_and_normalize_payload(payload: dict) -> dict:
    # Ensure required top-level keys
    required_keys = {"version", "input", "render", "style_preset", "annotations", "3d"}
    missing = [k for k in sorted(required_keys) if k not in payload]
    if missing:
        raise HTTPException(status_code=400, detail=f"Missing required keys: {', '.join(missing)}")

    if payload.get("version") != 1:
        raise HTTPException(status_code=400, detail="version must be 1")

    # Basic shape checks (lightweight)
    input_obj = payload.get("input", {})
    if input_obj.get("kind") not in ("smiles", "pdb", "name"):
        raise HTTPException(status_code=400, detail="input.kind must be 'smiles'|'pdb'|'name'")
    if not isinstance(input_obj.get("value"), str) or not input_obj.get("value"):
        raise HTTPException(status_code=400, detail="input.value must be non-empty string")
    if input_obj.get("conformer_method") not in ("etkdg", "none"):
        raise HTTPException(status_code=400, detail="input.conformer_method must be 'etkdg'|'none'")

    render_obj = payload.get("render", {})
    if not isinstance(render_obj.get("modes"), list) or not render_obj.get("modes"):
        raise HTTPException(status_code=400, detail="render.modes must be non-empty array")
    if not isinstance(render_obj.get("outputs"), list) or not render_obj.get("outputs"):
        raise HTTPException(status_code=400, detail="render.outputs must be non-empty array")
    for k in ("width", "height", "dpi"):
        v = render_obj.get(k)
        if not isinstance(v, (int, float)) or v <= 0:
            raise HTTPException(status_code=400, detail=f"render.{k} must be > 0")
    if not isinstance(render_obj.get("transparent"), bool):
        raise HTTPException(status_code=400, detail="render.transparent must be boolean")

    # 3d key must exist with exact spelling
    if "3d" not in payload or not isinstance(payload["3d"], dict):
        raise HTTPException(status_code=400, detail='Top-level key "3d" must be present and an object')

    return payload


def _normalize_for_hashing(payload: dict) -> dict:
    # Ensure the top-level key is exactly "3d". If model aliasing used, remap.
    obj = dict(payload)
    if "_3d" in obj and "3d" not in obj:
        obj["3d"] = obj.pop("_3d")
    return obj


def _start_background_render(spec: dict, spec_id: str) -> None:
    def _worker():
        # Update -> processing
        with _records_lock:
            rec = _records.get(spec_id)
            if not rec:
                return
            rec.status = "processing"

        try:
            _render_assets(spec, spec_id)
            urls = _asset_urls(spec_id)
            with _records_lock:
                rec = _records.get(spec_id)
                if rec:
                    rec.urls = urls or None
                    rec.status = "completed"
        except Exception as exc:  # pragma: no cover
            with _records_lock:
                rec = _records.get(spec_id)
                if rec:
                    rec.status = "failed"
                    rec.error = str(exc)

    t = threading.Thread(target=_worker, name=f"figure-render-{spec_id[:8]}", daemon=True)
    t.start()


@router.post("")
async def create_figure(spec: dict) -> dict:
    # Validate and normalize
    spec = _validate_and_normalize_payload(spec)
    normalized = _normalize_for_hashing(spec)

    # Compute spec_id deterministically
    sid = compute_spec_id(normalized)

    # Fast idempotent response if known
    with _records_lock:
        existing = _records.get(sid)
        if existing:
            response: Dict[str, Any] = {"spec_id": sid, "status": existing.status}
            if existing.urls is not None:
                response["urls"] = existing.urls
            return response

        # Create new record as queued and kick off background work
        _records[sid] = FigureRecord(spec_id=sid, status="queued", urls=None)

    _start_background_render(normalized, sid)

    return {"spec_id": sid, "status": "queued"}


@router.get("/{spec_id}")
async def get_figure(spec_id: str) -> dict:
    with _records_lock:
        rec = _records.get(spec_id)
        if not rec:
            raise HTTPException(status_code=404, detail="spec_id not found")
        response: Dict[str, Any] = {"spec_id": rec.spec_id, "status": rec.status}
        if rec.urls is not None:
            response["urls"] = rec.urls
        return response
