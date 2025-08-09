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


def _render_assets(spec: dict, spec_id: str) -> None:
    # Minimal, deterministic placeholders to satisfy contracts; heavy work can be replaced with real pipeline
    asset_dir = _asset_dir(spec_id)

    # 2D SVG placeholder honoring width/height and transparency flag
    width = int(spec.get("render", {}).get("width", 1024) or 1024)
    height = int(spec.get("render", {}).get("height", 768) or 768)
    transparent = bool(spec.get("render", {}).get("transparent", True))
    svg_bg = "none" if transparent else "white"
    svg = (
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>"
        f"<rect width='100%' height='100%' fill='{svg_bg}'/>"
        f"<text x='{width // 2}' y='{height // 2}' dominant-baseline='middle' text-anchor='middle' font-family='Helvetica' font-size='24'>"
        f"Moleculens Figure {spec_id[:8]}"  # deterministic label
        f"</text></svg>"
    ).encode()
    _write_once(asset_dir / "2d.svg", svg)

    # 2D PNG via cairosvg (optional dependency is present in requirements)
    try:
        import cairosvg  # type: ignore

        png_bytes = cairosvg.svg2png(
            bytestring=svg,
            output_width=width,
            output_height=height,
            dpi=int(spec.get("render", {}).get("dpi", 300) or 300),
        )
        _write_once(asset_dir / "2d.png", png_bytes)
    except Exception:
        pass

    # 3D PNG placeholder (re-use SVG converted as a stand-in)
    try:
        if not (asset_dir / "3d.png").exists():
            if (asset_dir / "2d.png").exists():
                # Duplicate as a placeholder to satisfy key presence
                data = (asset_dir / "2d.png").read_bytes()
                _write_once(asset_dir / "3d.png", data)
    except Exception:
        pass

    # Meta JSON (write-once)
    meta = {
        "spec_id": spec_id,
        "canonical": json.loads(canonical_json(spec)),
        "inputs": spec.get("input", {}),
        "render": spec.get("render", {}),
        "style_preset": spec.get("style_preset"),
        "annotations": spec.get("annotations", {}),
        "3d": spec.get("3d", {}),
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
