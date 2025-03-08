import os
import json
from agent_management.agents.pubchem_agent import PubChemAgent
from agent_management.llm_service import LLMService, LLMModelConfig, ProviderType

def test_save_compound_details():
    # Initialize the agent
    config = LLMModelConfig(
        provider=ProviderType.OPENAI,
        model_name="gpt-3.5-turbo",
        api_key=None
    )
    llm_service = LLMService(config)
    agent = PubChemAgent(llm_service)
    
    # Test CID
    cid = 4140276
    
    try:
        # Save compound details
        json_path = agent.save_compound_details_to_json(cid)
        
        # Verify the file exists
        assert os.path.exists(json_path), f"JSON file not found at {json_path}"
        
        # Read and print the contents
        with open(json_path, 'r') as f:
            data = json.load(f)
            print("\nCompound Details:")
            print("=" * 40)
            print(f"Name: {data['basic_info']['name']}")
            print(f"Formula: {data['basic_info']['molecular_formula']}")
            print(f"Weight: {data['basic_info']['molecular_weight']}")
            print(f"SMILES: {data['structure']['canonical_smiles']}")
            print("\nComposition:")
            print(f"Atoms: {data['composition']['atoms_count']}")
            print(f"Bonds: {data['composition']['bonds_count']}")
            print(f"Elements: {data['composition']['elements']}")
            print("\nProperties:")
            print(f"H-Bond Donors: {data['properties']['h_bond_donor_count']}")
            print(f"H-Bond Acceptors: {data['properties']['h_bond_acceptor_count']}")
            print(f"TPSA: {data['properties']['tpsa']}")
            print(f"XLogP: {data['properties']['xlogp']}")
            
            # Print full JSON structure
            print("\nFull JSON Structure:")
            print(json.dumps(data, indent=2))
            
    except Exception as e:
        print(f"Error: {str(e)}")
        raise

if __name__ == "__main__":
    test_save_compound_details() 