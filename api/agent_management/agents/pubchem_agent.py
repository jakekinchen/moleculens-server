"""
PubChem Agent - Handles interaction with PubChem database and molecular structure retrieval.
"""

import os
import tempfile
import json
import pubchempy as pcp
import requests
from typing import Optional, Dict, Any, NamedTuple
from agent_management.molecule_visualizer import MoleculeVisualizer
from agent_management.models import PubChemSearchResult, PubChemCompound
from agent_management.agents.script_agent import ScriptAgent
from agent_management.agents.pubchem_agent_helper import validate_and_convert_script
from agent_management.llm_service import (
    LLMService,
    LLMRequest,
    StructuredLLMRequest,
    LLMModelConfig,
    ProviderType
)
from rdkit import Chem
from rdkit.Chem import Fragments, Descriptors, AllChem

# Debug flag - set to False to disable debug logging
DEBUG_PUBCHEM = True

def _sdf_to_pdb_block(sdf_data: str) -> str:
        """
        Convert SDF data (string) to a single PDB block using RDKit in-memory.

        Returns an empty string if conversion fails.
        """
        mol = Chem.MolFromMolBlock(sdf_data, sanitize=True, removeHs=False)
        if mol is None:
            return ""

        if mol.GetNumConformers() == 0:
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())

        AllChem.MMFFOptimizeMolecule(mol)

        pdb_data = Chem.MolToPDBBlock(mol)
        return pdb_data if pdb_data else ""

def write_debug_file(filename: str, content: str) -> None:
    """Write debug content to a file if DEBUG_PUBCHEM is True."""
    if not DEBUG_PUBCHEM:
        return
        
    # Create debug directory if it doesn't exist
    debug_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "debug")
    os.makedirs(debug_dir, exist_ok=True)
    
    # Write content to file (overwriting if exists)
    filepath = os.path.join(debug_dir, filename)
    try:
        with open(filepath, 'w') as f:
            f.write(content)
        print(f"Debug file written: {filepath}")
    except Exception as e:
        print(f"Error writing debug file {filepath}: {e}")

class MoleculePackage(NamedTuple):
    """Container for molecule visualization package data"""
    js: str
    html: str
    title: str

class PubChemAgent:
    def __init__(self, llm_service: LLMService, use_element_labels: bool = True, 
                 convert_back_to_indices: bool = False, script_model: Optional[str] = None):
        self.llm_service = llm_service
        self.use_element_labels = use_element_labels  # Use element-based labels (C1, O1) by default
        self.convert_back_to_indices = convert_back_to_indices  # Convert back to numeric indices after script generation
        self.script_model = script_model  # Optional model override for script agent

    def interpret_user_query(self, user_input: str) -> str:
        """
        Use LLM to interpret user input into a molecule name or identifier.
        """
        request = LLMRequest(
            user_prompt=f"""Given this user input is related to molecular structures or trying to learn something about the microscopic structures of molecules or some topic related, give us the molecule that would best help the user learn what they are trying to learn about, return 'N/A' if you can't determine a valid molecule.
            User input: '{user_input}'
            Only respond with the molecule name or 'N/A', no other text.""",
            system_prompt="You are a chemistry professor that helps identify molecule names from user queries. For instance, if the user asks about wanting to learn about chiral carbon centers then you should return '2-Butanol' or 'CC(O)CC' and if the user says 'Teach me about transition metal complexes with octahedral geometry' then you should return 'titanocene dichloride' or Hexamminecobalt(III) chloride"
        )
        response = self.llm_service.generate(request)
        return response.content.strip()

    def fetch_sdf_for_cid(self, cid: int) -> Optional[str]:
        """
        Fetches the SDF string for a given PubChem CID using a direct REST call.
        Returns SDF as text, or None if something fails.
        """
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF"
        resp = requests.get(url, timeout=30)
        if resp.status_code == 200:
            return resp.text
        return None

    def get_molecule_sdfs(self, user_input: str) -> PubChemSearchResult:
        """
        1. Uses LLM to interpret the user input into a molecule name or ID.
        2. Queries PubChem for up to 5 matches.
        3. Fetches each match's SDF.
        4. Returns a structured result with the data.
        """
        # Step 1: Interpret the user's input via LLM
        molecule_name = self.interpret_user_query(user_input)
        secondary_molecule_name = self.interpret_user_query(user_input+" but not "+molecule_name + "so try to find a similar molecule that would be relevant to the user's query")
        if molecule_name == "N/A":
            # If the LLM can't parse or guess a name, return an empty structure
            return PubChemSearchResult(
                query=user_input,
                interpreted_query=molecule_name,
                results=[]
            )

        # Step 2: Query PubChem for the interpreted name
        compounds = pcp.get_compounds(molecule_name, 'name')
        if not compounds:
            compounds = pcp.get_compounds(secondary_molecule_name, 'name')
            if not compounds:
                return PubChemSearchResult(
                query=user_input,
                interpreted_query=molecule_name,
                results=[]
            )

        # Step 3: Take up to 5 matching CIDs
        limited_compounds = compounds[:5]

        # Step 4: For each CID, gather metadata + SDF
        results = []
        for cmpd in limited_compounds:
            sdf_str = self.fetch_sdf_for_cid(cmpd.cid)
            
            # Transform atoms and bonds into expected format
            atom_numbers = [atom.number for atom in cmpd.atoms] if cmpd.atoms else None
            bond_tuples = [(bond.aid1, bond.aid2, bond.order) for bond in cmpd.bonds] if cmpd.bonds else None
            
            results.append(PubChemCompound(
                name=molecule_name,
                cid=cmpd.cid,
                molecular_formula=cmpd.molecular_formula,
                molecular_weight=cmpd.molecular_weight,
                iupac_name=cmpd.iupac_name,
                sdf=sdf_str,
                canonical_smiles=cmpd.canonical_smiles,
                isomeric_smiles=cmpd.isomeric_smiles,
                elements=cmpd.elements,
                atoms=atom_numbers,
                bonds=bond_tuples,
                charge=cmpd.charge,
                synonyms=cmpd.synonyms
            ))

        # Step 5: Return structured result
        return PubChemSearchResult(
            query=user_input,
            interpreted_query=molecule_name,
            results=results
        )

    def get_compound_details(self, cid: int) -> Optional[PubChemCompound]:
        """
        Get detailed information for a specific compound by CID.
        """
        try:
            compound = pcp.Compound.from_cid(cid)
            sdf_str = self.fetch_sdf_for_cid(cid)
            
            return PubChemCompound(
                name=compound.iupac_name,  # Use IUPAC name as the primary name
                cid=compound.cid,
                molecular_formula=compound.molecular_formula,
                molecular_weight=compound.molecular_weight,
                iupac_name=compound.iupac_name,
                sdf=sdf_str,
                canonical_smiles=compound.canonical_smiles,
                isomeric_smiles=compound.isomeric_smiles,
                elements=compound.elements,
                atoms=compound.atoms,
                bonds=compound.bonds,
                charge=compound.charge,
                synonyms=compound.synonyms

                
                
            )
        except Exception as e:
            print(f"Error fetching compound details for CID {cid}: {str(e)}")
            return None
            
    def get_molecule_package(self, user_query: str) -> MoleculePackage:
        """
        Generate a complete molecule visualization package from a user query.
        
        Args:
            user_query (str): The user's natural language query about a molecule
            
        Returns:
            MoleculePackage: A tuple containing the JS code, HTML visualization, and title
            
        Raises:
            ValueError: If no valid molecules were found for the query
        """
        print(f"[DEBUG] Processing user query: {user_query}")
        
        try:
            # Step 1: Get molecule data from PubChem based on the query
            search_result = self.get_molecule_sdfs(user_query)
            print(f"[DEBUG] Search result interpreted query: {search_result.interpreted_query}")
            print(f"[DEBUG] Number of results: {len(search_result.results)}")
            
            # Step 2: Check if we have any results
            if not search_result.results:
                raise ValueError(f"No molecules found for query: '{user_query}' (interpreted as '{search_result.interpreted_query}')")
            
            # Step 3: Take the first result as our target molecule
            compound = search_result.results[0]
            print(f"[DEBUG] Selected compound: {compound.name} (CID: {compound.cid})")
            print(f"[DEBUG] IUPAC name: {compound.iupac_name}")
            
            if DEBUG_PUBCHEM:
                write_debug_file('pubchem_sdf.txt', compound.sdf or '')
            
            if compound.sdf is None:
                raise ValueError(f"No SDF data available for compound {compound.name} (CID: {compound.cid})")
            
            try:
                # Generate a display title for the molecule
                display_title = compound.name if compound.name else compound.iupac_name if compound.iupac_name else "Molecule"
                print(f"[DEBUG] Using display title: {display_title}")

                if compound.isomeric_smiles:
                    mol = Chem.MolFromSmiles(compound.isomeric_smiles)
                    if mol is None:
                        raise ValueError(f"Failed to parse SMILES for compound {compound.name} (CID: {compound.cid})")
                    
                    smarts_pattern = Chem.MolToSmarts(mol)

                    # NEW: Identify functional groups using RDKit fragment functions
                    functional_groups = {}
                    try:
                        functional_groups['aliphatic_OH_count'] = Fragments.fr_Al_OH(mol)
                        functional_groups['aromatic_OH_count'] = Fragments.fr_Ar_OH(mol)
                        functional_groups['halogen_count'] = Fragments.fr_halogen(mol)
                        functional_groups['amine_count'] = Fragments.fr_Al_N(mol)  # Changed from fr_Al_NH2 to fr_Al_N
                        print(f"[DEBUG] Functional groups: {functional_groups}")
                    except Exception as e:
                        print(f"[WARNING] Error counting functional groups: {str(e)}")
                        functional_groups = {
                            'aliphatic_OH_count': 0,
                            'aromatic_OH_count': 0,
                            'halogen_count': 0,
                            'amine_count': 0
                        }

                    # Create a list of functional groups to highlight
                    highlight_groups = []
                    if functional_groups['aliphatic_OH_count'] > 0:
                        highlight_groups.append('aliphatic_OH')
                    if functional_groups['aromatic_OH_count'] > 0:
                        highlight_groups.append('aromatic_OH')
                    if functional_groups['halogen_count'] > 0:
                        highlight_groups.append('halogen')
                    if functional_groups['amine_count'] > 0:
                        highlight_groups.append('amine')
                    print(f"[DEBUG] Highlight groups: {highlight_groups}")


                molecule_data = {
                    'name': compound.name,
                    'cid': compound.cid,
                    'smiles': compound.isomeric_smiles,
                    'highlight_groups': highlight_groups,
                    'smarts_pattern': smarts_pattern,
                    'functional_groups': functional_groups,
                    'iupac_name': compound.iupac_name,
                    'molecular_formula': compound.molecular_formula,
                    'molecular_weight': compound.molecular_weight,
                    'elements': compound.elements,
                    'atoms': compound.atoms,
                    'bonds': compound.bonds,
                    'charge': compound.charge,
                    'synonyms': compound.synonyms
                }

                print(f"[DEBUG] Molecule data: {molecule_data}")
                    

                # Send the molecule info and the user query to the script agent to get a script
                # Use script_model override if provided, otherwise use the same LLM service as PubChem agent
                if self.script_model:
                    from agent_management.agent_factory import AgentFactory
                    script_agent = AgentFactory.create_script_agent(self.script_model)
                else:
                    script_agent = ScriptAgent(self.llm_service)
                script = script_agent.generate_script_from_molecule(compound.name, user_query, molecule_data)
                
                # First convert to element-based labels for LLM understanding (if enabled)
                script = validate_and_convert_script(
                    script=script,
                    molecule_data=molecule_data,
                    use_element_labels=self.use_element_labels
                )
                
                # Convert back to numeric indices if configured that way
                if self.convert_back_to_indices:
                    script = validate_and_convert_script(
                        script=script,
                        molecule_data=molecule_data,
                        use_element_labels=False,
                        convert_back_to_indices=True
                    )
                
                print(f"[DEBUG] Generated script: {script}")

                pdb_data = _sdf_to_pdb_block(compound.sdf)
                
                # Generate the HTML content
                print("[DEBUG] Generating HTML content...")
                html_content = MoleculeVisualizer.generate_html_viewer_from_pdb(pdb_data, display_title)
                
                if DEBUG_PUBCHEM:
                    write_debug_file('pubchem_html.html', html_content)

                print(f"[DEBUG] Creating interactive visualization...")
                # Create a script data structure for the interactive visualization
                # The script variable appears to be an object with its own 'content' property,
                # but the JavaScript code expects scriptData.content to be an array directly
                if isinstance(script, dict) and 'content' in script:
                    # If script is a dict with a 'content' property, use that directly
                    script_data = {
                        "title": display_title,
                        "content": script['content']
                    }
                else:
                    # Otherwise, use the script as the content
                    script_data = {
                        "title": display_title,
                        "content": script
                    }
                
                print(f"[DEBUG] Script data structure: {type(script_data['content'])}")
                
                # Generate interactive HTML with both PDB data and script data
                interactive_html = MoleculeVisualizer.generate_interactive_html(
                    pdb_data=pdb_data, 
                    title=display_title,
                    script_data=script_data
                )
                
                if DEBUG_PUBCHEM:
                    write_debug_file('pubchem_interactive.html', interactive_html)
                
                # Generate minimal JS for embedding
                print("[DEBUG] Generating JS content...")
                js_content = MoleculeVisualizer.generate_js_code_from_pdb(pdb_data, display_title)
                
                if DEBUG_PUBCHEM:
                    write_debug_file('pubchem_js.js', js_content)
                
                # Create and log the package
                package = MoleculePackage(
                    js=js_content,
                    html=interactive_html,
                    title=display_title
                )
                
                if DEBUG_PUBCHEM:
                    debug_package = {
                        'title': package.title,
                        'js_length': len(package.js),
                        'html_length': len(package.html),
                        'sdf_length': len(compound.sdf)
                    }
                    write_debug_file('pubchem_package.json', json.dumps(debug_package, indent=2))
                
                print("[DEBUG] Successfully created molecule package")
                return package
            
            except Exception as e:
                print(f"[ERROR] Failed to generate visualization: {str(e)}")
                print(f"[ERROR] Error type: {type(e).__name__}")
                import traceback
                print(f"[ERROR] Traceback:\n{traceback.format_exc()}")
                raise ValueError(f"Failed to generate visualization: {str(e)}")
        
        except Exception as e:
            print(f"[ERROR] Error in get_molecule_package: {str(e)}")
            print(f"[ERROR] Error type: {type(e).__name__}")
            import traceback
            print(f"[ERROR] Traceback:\n{traceback.format_exc()}")
            raise