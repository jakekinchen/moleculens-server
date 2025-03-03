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
from agent_management.llm_service import (
    LLMService,
    LLMRequest,
    StructuredLLMRequest,
    LLMModelConfig,
    ProviderType
)

# Debug flag - set to False to disable debug logging
DEBUG_PUBCHEM = True

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
    def __init__(self, llm_service: LLMService):
        self.llm_service = llm_service

    def interpret_user_query(self, user_input: str) -> str:
        """
        Use LLM to interpret user input into a molecule name or identifier.
        """
        request = LLMRequest(
            user_prompt=f"""Given this user input is related to molecular structures or trying to learn something about the microscopic structures of molecules or some topic related, give us the molecule that would best help the user learn what they are trying to learn about, return 'N/A' if you can't determine a valid molecule.
            User input: '{user_input}'
            Only respond with the molecule name or 'N/A', no other text.""",
            system_prompt="You are a chemistry professor that helps identify molecule names from user queries. For instance, if the user asks about wanting to learn about chiral carbon centers then you should return '2-Butanol' or 'CC(O)CC'"
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
            results.append(PubChemCompound(
                name=molecule_name,
                cid=cmpd.cid,
                molecular_formula=cmpd.molecular_formula,
                molecular_weight=cmpd.molecular_weight,
                iupac_name=cmpd.iupac_name,
                sdf=sdf_str
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
                sdf=sdf_str
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
                display_title = compound.iupac_name if compound.iupac_name else compound.name
                print(f"[DEBUG] Using display title: {display_title}")
                
                # Generate the HTML content
                print("[DEBUG] Generating HTML content...")
                html_content = MoleculeVisualizer.generate_html_viewer_from_sdf(compound.sdf, display_title)
                
                if DEBUG_PUBCHEM:
                    write_debug_file('pubchem_html.html', html_content)
                
                # Generate minimal JS for embedding
                print("[DEBUG] Generating JS content...")
                js_content = MoleculeVisualizer.generate_js_code_from_sdf(compound.sdf, display_title)
                
                if DEBUG_PUBCHEM:
                    write_debug_file('pubchem_js.js', js_content)
                
                # Create and log the package
                package = MoleculePackage(
                    js=js_content,
                    html=html_content,
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