from typing import Dict, List, Optional

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

from api.pymol.services.pubchem import PubChemSearchService
from api.pymol.services.rdkit_utils import identify_functional_groups, smiles_to_mol

router = APIRouter(prefix="/chem", tags=["Chem"])


class FunctionalGroupRequest(BaseModel):
    cid: Optional[int] = Field(default=None, description="PubChem CID")
    smiles: Optional[str] = Field(default=None, description="SMILES string")
    include_matches: bool = Field(default=False, description="Include atom indices for each match")


class FunctionalGroupMatch(BaseModel):
    name: str
    atom_indices: List[int]


class FunctionalGroupResponse(BaseModel):
    groups: Dict[str, int]
    matches: Optional[List[FunctionalGroupMatch]] = None


@router.post("/functional-groups", response_model=FunctionalGroupResponse)
def functional_groups(payload: FunctionalGroupRequest) -> FunctionalGroupResponse:
    # Validate input: exactly one of cid or smiles must be provided
    if (payload.cid is None and not payload.smiles) or (payload.cid is not None and payload.smiles):
        raise HTTPException(status_code=400, detail="Provide exactly one of 'cid' or 'smiles'.")

    smiles: Optional[str] = payload.smiles

    if payload.cid is not None:
        service = PubChemSearchService()
        compound = service.get_compound_details(payload.cid)
        if compound is None:
            raise HTTPException(status_code=404, detail=f"No compound found for CID {payload.cid}.")
        # Try isomeric_smiles first, then canonical_smiles
        smiles = getattr(compound, "isomeric_smiles", None) or getattr(compound, "canonical_smiles", None)
        if not smiles:
            raise HTTPException(status_code=404, detail=f"No SMILES available for CID {payload.cid}.")

    assert smiles is not None  # for type checkers
    mol = smiles_to_mol(smiles)
    if mol is None:
        raise HTTPException(status_code=400, detail="Invalid SMILES string.")

    group_counts, match_list = identify_functional_groups(mol, include_matches=payload.include_matches)

    response = FunctionalGroupResponse(groups=group_counts)
    if payload.include_matches:
        response.matches = [FunctionalGroupMatch(name=name, atom_indices=list(indices)) for name, indices in match_list]
    return response
