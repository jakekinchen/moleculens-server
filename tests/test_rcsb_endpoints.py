from fastapi.testclient import TestClient
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import types
chem_mod = types.ModuleType('rdkit.Chem')
chem_mod.Fragments = object
chem_mod.Descriptors = object
chem_mod.AllChem = object
rdkit_mod = types.ModuleType('rdkit')
rdkit_mod.Chem = chem_mod
sys.modules['rdkit'] = rdkit_mod
sys.modules['rdkit.Chem'] = chem_mod
from api.main import app

client = TestClient(app)


def test_fetch_structure_pdb():
    resp = client.post('/rcsb/fetch-structure/', json={'identifier': '1STP', 'format': 'pdb'})
    assert resp.status_code == 200
    data = resp.json()
    assert 'ATOM' in data['data']
