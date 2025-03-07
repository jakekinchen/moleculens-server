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
import logging
from agent_management.debug_utils import DEBUG_PUBCHEM, write_debug_file
import datetime

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

class MoleculePackage(NamedTuple):
    """Container for molecule visualization package data"""
    pdb_data: str
    html: str
    title: str

class PubChemAgent:
    """Agent for retrieving and processing molecular data from PubChem"""
    
    def __init__(self, llm_service: LLMService, use_element_labels: bool = True, 
                 convert_back_to_indices: bool = False, script_model: Optional[str] = None):
        """
        Initialize the PubChem agent
        
        Args:
            llm_service: LLM service for generating interpretations and scripts
            use_element_labels: Whether to use element-based labels (C1, O1) instead of indices
            convert_back_to_indices: Whether to convert element-labels back to numeric indices
            script_model: Optional model override for script agent
        """
        self.llm_service = llm_service
        self.use_element_labels = use_element_labels
        self.convert_back_to_indices = convert_back_to_indices  # Convert back to numeric indices after script generation
        self.script_model = script_model  # Optional model override for script agent
        self.logger = logging.getLogger(__name__)

    def interpret_user_query(self, user_input: str) -> str:
        """
        Use LLM to interpret user input into a molecule name or identifier.
        Handles provider failures gracefully with fallback behavior.
        
        Args:
            user_input: The user's query about a molecule
            
        Returns:
            A molecule name or identifier, or the original input if interpretation fails
            
        Raises:
            Exception: If all interpretation attempts fail and the input cannot be used directly
        """
        try:
            request = LLMRequest(
                user_prompt=f"""Given this user input is related to molecular structures or trying to learn something about the microscopic structures of molecules or some topic related, give us the molecule that would best help the user learn what they are trying to learn about, return 'N/A' if you can't determine a valid molecule.
                User input: '{user_input}'
                Only respond with the molecule name or 'N/A', no other text.""",
                system_prompt="You are a chemistry professor that helps identify molecule names from user queries. For instance, if the user asks about wanting to learn about chiral carbon centers then you should return '2-Butanol' or 'CC(O)CC' and if the user says 'Teach me about transition metal complexes with octahedral geometry' then you should return 'titanocene dichloride' or Hexamminecobalt(III) chloride"
            )
            response = self.llm_service.generate(request)
            return response.content.strip()
        except Exception as e:
            self.logger.warning(f"Failed to interpret user query with LLM: {str(e)}")
            
            # Fallback 1: Try to use the input directly if it looks like a molecule name
            if len(user_input.split()) <= 3 and not any(char in user_input for char in "?!.,:;"):
                self.logger.info(f"Using user input directly as molecule name: '{user_input}'")
                return user_input
                
            # Fallback 2: Use a default molecule if the input is complex
            self.logger.info(f"Using default molecule 'water' for complex query: '{user_input}'")
            return "water"

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

    def get_molecule_sdfs(self, user_input: str, skip_llm: bool = False) -> PubChemSearchResult:
        """
        1. Uses LLM to interpret the user input into a molecule name or ID.
        2. Queries PubChem for up to 5 matches.
        3. Fetches each match's SDF.
        4. Returns a structured result with the data.
        """
        if skip_llm:
            # If skip_llm is True, directly use the input as the molecule name
            molecule_name = user_input
        else:
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
        Get a complete molecule visualization package from a user query.
        
        Args:
            user_query: The user's query about a molecule
            
        Returns:
            MoleculePackage: A package containing PDB data, HTML, and title for visualization
            
        Raises:
            ValueError: If no molecules are found or other errors occur
        """
        # Write a debug file to confirm this method is being called
        if DEBUG_PUBCHEM:
            debug_info = {
                'method': 'get_molecule_package',
                'user_query': user_query,
                'timestamp': datetime.datetime.now().isoformat()
            }
            write_debug_file('pubchem_package_start.json', json.dumps(debug_info, indent=2))
            
        self.logger.info(f"[DEBUG] Getting molecule package for query: {user_query}")
        
        try:
            # First, get the molecule SDFs
            search_result = self.get_molecule_sdfs(user_query)
            
            if not search_result.results:
                raise ValueError(f"No molecules found for query: {user_query}, search result: {search_result}")
                
            # Use the first compound
            compound = search_result.results[0]

            if compound.sdf is None:
                raise ValueError(f"No SDF data available for compound {compound.name} (CID: {compound.cid})")

            pdb_data = _sdf_to_pdb_block(compound.sdf)
            
            # Get a display title (use name if available, otherwise ID)
            display_title = compound.name if compound.name else f"CID {compound.cid}"
            self.logger.info(f"[DEBUG] Using molecule: {display_title}")

            if DEBUG_PUBCHEM:
                write_debug_file('pubchem_sdf.txt', compound.sdf or '')

            try:
                # Generate a display title for the molecule
                display_title = compound.name if compound.name else compound.iupac_name if compound.iupac_name else "Molecule"
                self.logger.info(f"Using display title: {display_title}")

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
                        self.logger.info(f"[DEBUG] Functional groups: {functional_groups}")
                    except Exception as e:
                        self.logger.warning(f"[WARNING] Error counting functional groups: {str(e)}")
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
                    self.logger.info(f"[DEBUG] Highlight groups: {highlight_groups}")


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

                self.logger.info(f"[DEBUG] Molecule data: {molecule_data}")
                    

                # Send the molecule info and the user query to the script agent to get a script
                # Use script_model override if provided, otherwise use the same LLM service as PubChem agent
                if self.script_model:
                    from agent_management.agent_factory import AgentFactory
                    script_agent = AgentFactory.create_script_agent(self.script_model)
                else:
                    script_agent = ScriptAgent(self.llm_service)

                # First convert molecule data to use elemental indices for better LLM understanding
                from agent_management.agents.pubchem_agent_helper import generate_atom_label_mapping
                atom_label_mapping = generate_atom_label_mapping(molecule_data)
                
                # Update molecule_data with element labels for the script agent
                molecule_data_with_labels = molecule_data.copy()
                molecule_data_with_labels['atom_labels'] = atom_label_mapping
                
                # Generate script using molecule data with elemental indices
                script = script_agent.generate_script_from_molecule(compound.name, user_query, molecule_data_with_labels)
                
                self.logger.info(f"[DEBUG] Generated raw script: {script}")
                
                # We ALWAYS convert elemental indices (C1, H1, etc.) back to numeric indices (0, 1, 2)
                # This is required for PDB visualization regardless of the use_element_labels setting
                script = validate_and_convert_script(
                    script=script,
                    molecule_data=molecule_data,
                    use_element_labels=False,
                    convert_back_to_indices=True
                )
                    
                self.logger.info(f"[DEBUG] Processed script with numeric indices: {script}")
                
                # Explicitly verify that all atom references are numeric indices, not element labels
                if DEBUG_PUBCHEM:
                    try:
                        all_numeric = True
                        for time_point in script['content']:
                            for atom_ref in time_point['atoms']:
                                # Check if any atom reference contains alphabet characters (would indicate element labels)
                                if any(c.isalpha() for c in str(atom_ref)):
                                    self.logger.error(f"[ERROR] Found non-numeric atom reference: {atom_ref}")
                                    all_numeric = False
                        
                        if all_numeric:
                            self.logger.info("[DEBUG] Verification passed: All atom references are numeric")
                        else:
                            self.logger.error("[ERROR] Verification failed: Some atom references are not numeric")
                    except Exception as e:
                        self.logger.error(f"[ERROR] Error during atom reference verification: {str(e)}")
                
                self.logger.info(f"[DEBUG] Generated script: {script}")

                
                
                # Generate the HTML content
                self.logger.info("[DEBUG] Generating HTML content...")
                html_content = MoleculeVisualizer.generate_html_viewer_from_pdb(pdb_data, display_title)
                
                if DEBUG_PUBCHEM:
                    write_debug_file('pubchem_html.html', html_content)

                self.logger.info(f"[DEBUG] Creating interactive visualization...")
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
                
                self.logger.info(f"[DEBUG] Script data structure: {type(script_data['content'])}")
                
                # Generate interactive HTML with both PDB data and script data
                interactive_html = MoleculeVisualizer.generate_interactive_html(
                    pdb_data=pdb_data, 
                    title=display_title,
                    script_data=script_data
                )
                
                if DEBUG_PUBCHEM:
                    write_debug_file('pubchem_interactive.html', interactive_html)
                
                # Generate minimal JS for embedding
                self.logger.info("[DEBUG] Generating JS content...")
                js_content = MoleculeVisualizer.generate_js_code_from_pdb(pdb_data, display_title)
                
                if DEBUG_PUBCHEM:
                    write_debug_file('pubchem_js.js', js_content)

                # Create the package
                self.logger.info("[DEBUG] Creating molecule package")
                package = MoleculePackage(
                    pdb_data=pdb_data,
                    html=interactive_html,
                    title=display_title
                )
                
                if DEBUG_PUBCHEM:
                    try:
                        debug_package = {
                            'title': package.title,
                            'pdb_data_length': len(package.pdb_data),
                            'html_length': len(package.html),
                            'sdf_length': len(compound.sdf),
                            'timestamp': datetime.datetime.now().isoformat()
                        }
                        write_debug_file('pubchem_package.json', json.dumps(debug_package, indent=2))
                        self.logger.info("[DEBUG] Package debug info written to: pubchem_package.json")
                    except Exception as e:
                        error_msg = f"Error writing package debug file: {str(e)}"
                        self.logger.error(f"[ERROR] {error_msg}")
                        write_debug_file('pubchem_package_write_error.txt', error_msg)
                
                self.logger.info("[DEBUG] Successfully created molecule package")
                return package
            
            except Exception as e:
                self.logger.error(f"[ERROR] Failed to generate visualization: {str(e)}")
                self.logger.error(f"[ERROR] Error type: {type(e).__name__}")
                import traceback
                self.logger.error(f"[ERROR] Traceback:\n{traceback.format_exc()}")
                raise ValueError(f"Failed to generate visualization: {str(e)}")
        
        except Exception as e:
            self.logger.error(f"[ERROR] Error in get_molecule_package: {str(e)}")
            self.logger.error(f"[ERROR] Error type: {type(e).__name__}")
            import traceback
            self.logger.error(f"[ERROR] Traceback:\n{traceback.format_exc()}")
            raise