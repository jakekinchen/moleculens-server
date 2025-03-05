"""
Comprehensive test script for the Groq provider with all response models
"""

import os
import sys
from typing import List, Dict, Optional, Any, Union
from pydantic import BaseModel, Field
from dotenv import load_dotenv

# Add the parent directory to the path so we can import the modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api.agent_management.llm_service import LLMService, StructuredLLMRequest, LLMModelConfig, ProviderType
from api.agent_management.models import (
    Vector3, Material, Geometry, Mesh, ThreeGroup, BooleanResponse, SceneObject
)

# Load environment variables
load_dotenv()

# Define additional test models
class Recipe(BaseModel):
    name: str
    ingredients: List[str]
    steps: List[str]
    prep_time_minutes: Optional[int] = None
    cook_time_minutes: Optional[int] = None

class MovieReview(BaseModel):
    title: str
    rating: float = Field(ge=0.0, le=10.0)
    review_text: str
    pros: List[str]
    cons: List[str]

def test_groq_all_models():
    """Test all response models with the Groq provider"""
    print("Testing All Response Models with Groq Provider...")
    
    # Check if GROQ_API_KEY is set
    if not os.getenv("GROQ_API_KEY"):
        print("Error: GROQ_API_KEY environment variable is not set")
        return False
    
    # Check if groq package is installed
    try:
        import groq
        print("Groq package is installed")
    except ImportError:
        print("Error: groq package is not installed")
        print("Install with: pip install groq")
        return False
    
    # Create an LLM service with Groq provider
    llm_service = LLMService(
        config=LLMModelConfig(
            provider=ProviderType.GROQ,
            model_name="llama3-70b-8192"
        )
    )
    
    # Test models
    test_results = {}
    
    # 1. Test Recipe model
    print("\n1. Testing Recipe model...")
    try:
        request = StructuredLLMRequest(
            user_prompt="Give me a recipe for a simple chocolate chip cookie.",
            system_prompt="You are a helpful cooking assistant that provides recipes in a structured format.",
            response_model=Recipe
        )
        
        recipe = llm_service.generate_structured(request)
        print(f"Recipe: {recipe.name}")
        print(f"Ingredients: {len(recipe.ingredients)} ingredients")
        print(f"Steps: {len(recipe.steps)} steps")
        test_results["Recipe"] = "✅ Passed"
    except Exception as e:
        print(f"Error: {str(e)}")
        test_results["Recipe"] = f"❌ Failed: {str(e)}"
    
    # 2. Test BooleanResponse model
    print("\n2. Testing BooleanResponse model...")
    try:
        request = StructuredLLMRequest(
            user_prompt="Is the following statement true? 'Water boils at 100 degrees Celsius at sea level.'",
            system_prompt="You are a fact-checking assistant. Determine if statements are true or false.",
            response_model=BooleanResponse
        )
        
        bool_response = llm_service.generate_structured(request)
        print(f"Is true: {bool_response.is_true}")
        print(f"Confidence: {bool_response.confidence}")
        print(f"Reasoning: {bool_response.reasoning[:100]}...")
        test_results["BooleanResponse"] = "✅ Passed"
    except Exception as e:
        print(f"Error: {str(e)}")
        test_results["BooleanResponse"] = f"❌ Failed: {str(e)}"
    
    # 3. Test MovieReview model
    print("\n3. Testing MovieReview model...")
    try:
        request = StructuredLLMRequest(
            user_prompt="Write a review for the movie 'The Matrix'.",
            system_prompt="You are a movie critic who writes structured reviews.",
            response_model=MovieReview
        )
        
        movie_review = llm_service.generate_structured(request)
        print(f"Movie: {movie_review.title}")
        print(f"Rating: {movie_review.rating}/10")
        print(f"Review: {movie_review.review_text[:100]}...")
        print(f"Pros: {len(movie_review.pros)} points")
        print(f"Cons: {len(movie_review.cons)} points")
        test_results["MovieReview"] = "✅ Passed"
    except Exception as e:
        print(f"Error: {str(e)}")
        test_results["MovieReview"] = f"❌ Failed: {str(e)}"
    
    # 4. Test Vector3 model
    print("\n4. Testing Vector3 model...")
    try:
        request = StructuredLLMRequest(
            user_prompt="Generate 3D coordinates for a point in space.",
            system_prompt="You are a 3D modeling assistant. Generate coordinates in 3D space.",
            response_model=Vector3
        )
        
        vector = llm_service.generate_structured(request)
        print(f"Vector3: ({vector.x}, {vector.y}, {vector.z})")
        test_results["Vector3"] = "✅ Passed"
    except Exception as e:
        print(f"Error: {str(e)}")
        test_results["Vector3"] = f"❌ Failed: {str(e)}"
    
    # 5. Test ThreeGroup model (complex nested model)
    print("\n5. Testing ThreeGroup model...")
    try:
        request = StructuredLLMRequest(
            user_prompt="Create a 3D model of a simple molecule with two atoms connected by a bond.",
            system_prompt="You are a 3D modeling expert. Generate precise Three.js structures.",
            response_model=ThreeGroup,
            temperature=0.2  # Low temperature for more deterministic output
        )
        
        three_group = llm_service.generate_structured(request)
        print(f"Group name: {three_group.name}")
        print(f"Number of meshes: {len(three_group.children)}")
        print(f"Position: ({three_group.position.x}, {three_group.position.y}, {three_group.position.z})")
        test_results["ThreeGroup"] = "✅ Passed"
    except Exception as e:
        print(f"Error: {str(e)}")
        test_results["ThreeGroup"] = f"❌ Failed: {str(e)}"
    
    # 6. Test SceneObject model
    print("\n6. Testing SceneObject model...")
    try:
        request = StructuredLLMRequest(
            user_prompt="Create a scene object representing a water molecule.",
            system_prompt="You are a scientific visualization expert. Create detailed scene objects.",
            response_model=SceneObject
        )
        
        scene_object = llm_service.generate_structured(request)
        print(f"Object name: {scene_object.name}")
        print(f"Description: {scene_object.description[:100]}...")
        print(f"Properties: {len(scene_object.properties)} properties")
        print(f"Appears at: {scene_object.appears_at}")
        test_results["SceneObject"] = "✅ Passed"
    except Exception as e:
        print(f"Error: {str(e)}")
        test_results["SceneObject"] = f"❌ Failed: {str(e)}"
    
    # Print summary of results
    print("\n" + "="*50)
    print("SUMMARY OF TEST RESULTS")
    print("="*50)
    for model, result in test_results.items():
        print(f"{model}: {result}")
    
    # Overall result
    if all("Passed" in result for result in test_results.values()):
        print("\n✅ All model tests passed with Groq provider")
        return True
    else:
        print("\n❌ Some model tests failed with Groq provider")
        return False

if __name__ == "__main__":
    test_groq_all_models() 