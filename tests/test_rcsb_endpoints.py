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

# Stub out the heavy `openai` dependency used in some modules
openai_mod = types.ModuleType('openai')
openai_mod.OpenAI = object
types_mod = types.ModuleType('openai.types')
chat_mod = types.ModuleType('openai.types.chat')
chat_completion_mod = types.ModuleType('openai.types.chat.chat_completion')
setattr(types_mod, 'Completion', object)
setattr(chat_mod, 'ChatCompletion', object)
setattr(chat_mod, 'ChatCompletionMessage', object)
setattr(chat_mod, 'ChatCompletionMessageParam', dict)
setattr(chat_completion_mod, 'Choice', object)
sys.modules['openai'] = openai_mod
sys.modules['openai.types'] = types_mod
sys.modules['openai.types.chat'] = chat_mod
sys.modules['openai.types.chat.chat_completion'] = chat_completion_mod
from api.main import app

client = TestClient(app)


def test_fetch_structure_pdb():
    resp = client.post('/rcsb/fetch-structure/', json={'identifier': '1STP', 'format': 'pdb'})
    assert resp.status_code == 200
    data = resp.json()
    assert 'ATOM' in data['data']
