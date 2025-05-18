"""
Script Agent - Generates structured scene scripts with timecodes, descriptions, and captions.
This is separate from the animation agent that will handle actual animation generation.
"""

from agent_management.models import SceneScript
from typing import Dict, Any
import json
from agent_management.llm_service import (
    LLMService,
    StructuredLLMRequest,
    LLMModelConfig,
    ProviderType
)

class ScriptAgent:
    def __init__(self, llm_service: LLMService):
        self.llm_service = llm_service
    
    def generate_script(self, topic: str) -> SceneScript:
        """
        Generate a structured scene script for a scientific topic.
        
        Args:
            topic: The scientific topic to explain in the scene
            
        Returns:
            SceneScript: A Pydantic model instance of the structured script
        """
        # Create a request for script generation
        request = StructuredLLMRequest(
            user_prompt=f"""Create a detailed script for a 3D scientific visualization about: {topic}
            
The script should cover a single coherent 3D scene that effectively explains the topic visually.
Focus on creating a cohesive visual narrative with clear transition points.

Return a JSON object with:
1. A concise, descriptive title for the scene
2. A content array containing 5-8 key points in the scene, each with:
   - timecode: Time marker in MM:SS format (starting at 00:00, ending around 02:00)
   - description: Detailed explanation of what should be happening in the 3D scene at this point
   - caption: An educational text caption that would appear on screen (50-100 characters)
   
CAPTION REQUIREMENTS:
- Each caption must be self-contained and meaningful on its own
- Directly relate to what is currently visible in the scene
- Use clear, concise language that balances technical accuracy with accessibility
- Highlight the key scientific concept being demonstrated at that moment
- Build understanding progressively through the scene

Example of good captions:
- "Carbon forms tetrahedral bonds, creating a 3D pyramid structure"
- "Water molecules bend at 104.5°, giving them a unique polar shape"
- "Electron clouds overlap as covalent bonds form between atoms"

Make sure each time point builds logically on the previous ones to tell a complete story about the topic.
Ensure descriptions are specific enough to guide the creation of 3D visuals (objects, movements, transitions).
Keep the explanations scientific but accessible to a general audience.""",
            
            system_prompt="""You are an expert in scientific communication and 3D visualization. Your task is to create structured scripts for educational scientific scenes.

When creating scene scripts:
- Focus on clear, accurate scientific explanations
- Create logical progression of visual elements
- Design effective transitions between key concepts
- Balance visual appeal with educational value
- Use appropriate pacing (approx. 2 minutes total)
- Structure time points to flow naturally from introduction to conclusion

CAPTION WRITING GUIDELINES:
1. Educational Value:
   - Each caption should teach a specific scientific concept
   - Use precise scientific terminology when relevant
   - Connect visual elements to underlying principles

2. Accessibility:
   - Write at a high school science level
   - Define complex terms within the caption when needed
   - Use active voice and present tense
   - Avoid jargon unless necessary for understanding

3. Visual Connection:
   - Describe what is currently visible on screen
   - Reference specific visual elements being highlighted
   - Use spatial terms to guide attention
   - Connect visual changes to scientific concepts

4. Progressive Learning:
   - Build complexity gradually through the scene
   - Reference previously introduced concepts when relevant
   - Create clear connections between sequential points
   - End with synthesis of key concepts

Example Caption Progression:
00:00 - "Methane molecule (CH4) with central carbon atom"
00:15 - "Four hydrogen atoms arrange in a tetrahedral pattern"
00:30 - "109.5° bond angles maximize electron cloud separation"
00:45 - "This tetrahedral shape gives methane its unique properties"

Always return complete, properly formatted JSON objects matching the requested schema exactly.
""",
            llm_config=self.llm_service.config,
            response_model=SceneScript,
        )
        
        script = self.llm_service.generate_structured(request)
        return script
    
    def generate_script_from_molecule(self, molecule_name: str, user_query: str, molecule_data: Dict[str, Any]) -> SceneScript:
        """
        Generate a structured scene script for a molecule.
        
        Returns:
            SceneScript: A Pydantic model instance of the structured script
        """
        # Check if we have atom labels in the molecule data
        atom_labels_info = ""
        if 'atom_labels' in molecule_data:
            # Extract the atom labels mapping and format it for better LLM understanding
            atom_labels = molecule_data['atom_labels']
            atom_labels_formatted = "\n".join([f"{idx}: {label}" for idx, label in atom_labels.items()])
            atom_labels_info = f"""
ATOM LABELS MAPPING:
The molecule has the following atom labels that you should use in your script:
{atom_labels_formatted}

When referring to atoms in your script, use these atom labels (like C1, H2, O1) 
rather than numeric indices, as these labels are more chemically meaningful.
"""
        
        # Convert molecule_data to string representation safely
        # Remove 'atom_labels' from the JSON to avoid redundancy
        molecule_data_copy = molecule_data.copy()
        if 'atom_labels' in molecule_data_copy:
            del molecule_data_copy['atom_labels']
        molecule_json = json.dumps(molecule_data_copy, indent=2)

        system_prompt = f"""You are an expert in scientific communication and 3D visualization. Your task is to create structured, engaging scripts for educational scientific scenes.

Given:
- A molecule name: {molecule_name}
- A user query indicating the area of interest: {user_query}
- Supplemental data about the molecule: {molecule_json}
{atom_labels_info}

Your Objective:
Create a clear, informative script that explains the molecule's structure and properties, guided by the user's area of interest.

Script Guidelines:
1. Clearly introduce the molecule using its IUPAC name (if applicable) and connect its structural geometry to VSEPR theory.
2. Explain naturally and logically, transitioning smoothly between key points.
3. Align the script's complexity and length to the molecule's complexity:
    - Simple molecules (e.g., water, methane): 5-7 key points (~1 min)
    - Complex molecules (e.g., glucose, insulin): 9-12 key points (~2 min)
4. Use the atoms field to highlight the specific atoms in the molecule that are relevant to the current caption's focus.

IMPORTANT FORMATTING REQUIREMENTS:
1. The "atoms" field MUST contain an array of STRING values, not numbers
2. If atom labels are provided above, use those labels (like "C1", "H2", "O1") in the atoms array
3. If no atom labels are provided, use string indices like "0", "1", "2"
4. NEVER include integers directly in the atoms array, ALWAYS wrap them in quotes
5. The introduction caption should have an empty atoms array.

Output Format:
Return a JSON object structured precisely as follows (realistic example):

{{
    "title": "Benzene: Aromaticity and Its Chemical Significance",
    "content": [
        {{
            "timecode": "00:00",
            "atoms": [],
            "caption": "Benzene's hexagonal ring is key to its aromatic stability."
        }},
        {{
            "timecode": "00:08",
            "atoms": ["C1", "C2", "C3", "C4", "C5", "C6"],
            "caption": "This unique geometry affects benzene's chemical reactivity."
        }},
        {{
            "timecode": "00:15",
            "atoms": ["C1", "C2", "C3", "C4", "C5", "C6", "H1", "H2", "H3", "H4", "H5", "H6"],
            "caption": "Symmetrical hydrogen placement reduces molecular polarity."
        }},
        {{
            "timecode": "00:22",
            "atoms": ["C1", "C2", "C3", "C4", "C5", "C6", "H1", "H2", "H3", "H4", "H5", "H6"],
            "caption": "Low polarity explains benzene's limited solubility in water."
        }},
        {{
            "timecode": "00:30",
            "atoms": ["C1", "C2", "C3", "C4", "C5", "C6"],
            "caption": "Electron delocalization creates stable aromatic pi bonds."
        }},
        {{
            "timecode": "00:37",
            "atoms": ["C1", "C2", "C3", "C4", "C5", "C6"],
            "caption": "Stable pi bonds influence benzene's resistance to addition reactions."
        }},
        {{
            "timecode": "00:45",
            "atoms": ["C1", "C2", "C3", "C4", "C5", "C6"],
            "caption": "Planar resonance structure contributes to overall molecular stability."
        }},
        {{
            "timecode": "00:52",
            "atoms": ["C1", "C2", "C3", "C4", "C5", "C6"],
            "caption": "This stability makes benzene an important industrial chemical."
        }},
        {{
            "timecode": "01:00",
            "atoms": ["C1", "C2", "C3", "C4", "C5", "C6"],
            "caption": "Equal bond lengths result from resonance and electron delocalization."
        }},
        {{
            "timecode": "01:07",
            "atoms": ["C1", "C2", "C3", "C4", "C5", "C6"],
            "caption": "Uniform bonding impacts benzene's predictable chemical behavior."
        }},
        {{
            "timecode": "01:15",
            "atoms": ["C1", "C2", "C3", "C4", "C5", "C6"],
            "caption": "Aromatic stability significantly influences benzene's reaction pathways."
        }},
        {{
            "timecode": "01:22",
            "atoms": ["C1", "C2", "C3", "C4", "C5", "C6"],
            "caption": "Understanding benzene's aromaticity helps explain its industrial uses."
        }}
    ]
}}

Use the user query to guide the emphasis of your script content without explicitly referencing it. For example, if the query is "Explain the aromatic rings in benzene," focus on aromaticity, molecular geometry, and its significance, incorporating other molecular details only as contextually relevant.
"""

        request = StructuredLLMRequest(
            user_prompt="Generate a molecule visualization script that follows the provided guidelines.",
            system_prompt=system_prompt,
            llm_config=self.llm_service.config,
            response_model=SceneScript,
        )

        script = self.llm_service.generate_structured(request)
        return script
