"""FastAPI routes for deterministic nomenclature explanations."""

from fastapi import APIRouter, HTTPException

from moleculens.api.schemas import (
    NomenclatureComponent,
    NomenclatureRequest,
    NomenclatureResponse,
    NomenclatureTraceStep,
)
from moleculens.compute.nomenclature import (
    NomenclatureError,
    UnsupportedNomenclatureError,
    explain_nomenclature,
)
from moleculens.core import get_logger

logger = get_logger(__name__)

router = APIRouter(prefix="/v1/nomenclature", tags=["nomenclature"])


@router.post("/explain", response_model=NomenclatureResponse)
def explain(request: NomenclatureRequest) -> NomenclatureResponse:
    """Explain a generated name and map referenced atoms to viewer order when possible."""

    try:
        explanation = explain_nomenclature(
            request.smiles,
            molblock=request.molblock,
            verify_opsin=request.verify_opsin,
        )
        return NomenclatureResponse(
            displayName=explanation.display_name,
            verified=explanation.verified,
            opsinStatus=explanation.opsin_status,
            rulesHit=explanation.rules_hit,
            ruleHints=explanation.rule_hints,
            trace=[
                NomenclatureTraceStep(
                    phase=step.phase,
                    decision=step.decision,
                    reason=step.reason,
                    atoms=step.atoms,
                    engineAtoms=step.engine_atoms,
                    data=step.data,
                )
                for step in explanation.trace
            ],
            components=[
                NomenclatureComponent(phase=component.phase, text=component.text)
                for component in explanation.components
            ],
            atomLocants=explanation.atom_locants,
            mappingStatus=explanation.mapping_status,
            engineVersion=explanation.engine_version,
            rdkitVersion=explanation.rdkit_version,
        )
    except UnsupportedNomenclatureError as exc:
        raise HTTPException(
            status_code=422,
            detail={"code": exc.code, "message": str(exc)},
        ) from exc
    except NomenclatureError as exc:
        raise HTTPException(
            status_code=422,
            detail={"code": "naming_failed", "message": str(exc)},
        ) from exc
    except Exception as exc:
        logger.exception("Nomenclature explanation failed", error=str(exc))
        raise HTTPException(
            status_code=500,
            detail={"code": "nomenclature_error", "message": "Nomenclature engine failed"},
        ) from exc
