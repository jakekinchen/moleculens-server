import os
import sys
import json
from typing import List, Dict, Optional
from pydantic import BaseModel, Field
from dotenv import load_dotenv

# Add the parent directory to sys.path to import modules from api
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import the GroqProvider and LLMRequest
from agent_management.providers.groq_provider import GroqProvider
from agent_management.llm_service import LLMRequest, LLMModelConfig, ProviderType

# Load environment variables
load_dotenv()

# Define a simple Pydantic model for structured output
class UserInfo(BaseModel):
    name: str
    age: int
    occupation: Optional[str] = None
    skills: List[str] = []

class Recipe(BaseModel):
    title: str
    ingredients: List[str]
    instructions: List[str]
    prep_time_minutes: int
    cook_time_minutes: int
    difficulty: str

def test_groq_structured_simple():
    """
    Test structured output generation with Groq without using instructor.
    """
    print("\n==================================================")
    print("Testing Groq for Structured Output (Simple Approach)")
    print("==================================================")
    
    # Check if GROQ_API_KEY is set
    api_key = os.environ.get("GROQ_API_KEY")
    if not api_key:
        print("Error: GROQ_API_KEY environment variable is not set.")
        print("Please set the GROQ_API_KEY environment variable and try again.")
        return False
    
    # Check if groq package is installed
    try:
        import groq
        print("Groq package is installed")
    except ImportError:
        print("Error: groq package is not installed.")
        print("Please install the groq package using: pip install groq")
        return False
    
    # Create a GroqProvider instance
    provider = GroqProvider(api_key=api_key)
    
    # Test 1: Extract user information using a structured prompt
    print("\nTest 1: Extract user information")
    user_text = """
    John Doe is a 35-year-old software engineer who specializes in Python, JavaScript, and cloud technologies.
    He has been working in the tech industry for over 10 years.
    """
    
    user_prompt = f"""
    Extract the following information from the text below and format it as a valid JSON object with these fields:
    - name: The person's full name
    - age: The person's age as an integer
    - occupation: The person's job or profession
    - skills: A list of the person's skills or technologies they know
    
    Text: {user_text}
    
    JSON:
    """
    
    # Create a proper LLMRequest
    request = LLMRequest(
        user_prompt=user_prompt,
        llm_config=LLMModelConfig(
            provider=ProviderType.GROQ,
            model_name="llama3-70b-8192",  # Using Llama 3 70B model
            api_key=api_key
        ),
        max_tokens=500
    )
    
    try:
        response = provider.generate(request)
        print("Raw response:")
        print(response.content)
        
        # Try to parse the JSON from the response
        try:
            # Find JSON in the response (it might be surrounded by markdown code blocks or other text)
            json_str = response.content
            if "```json" in response.content:
                json_str = response.content.split("```json")[1].split("```")[0].strip()
            elif "```" in response.content:
                json_str = response.content.split("```")[1].split("```")[0].strip()
            
            user_data = json.loads(json_str)
            user_info = UserInfo(**user_data)
            print("\nParsed UserInfo:")
            print(f"Name: {user_info.name}")
            print(f"Age: {user_info.age}")
            print(f"Occupation: {user_info.occupation}")
            print(f"Skills: {', '.join(user_info.skills)}")
            print("Test 1 passed!")
        except Exception as e:
            print(f"Error parsing JSON: {e}")
            print("Test 1 failed!")
            return False
    except Exception as e:
        print(f"Error generating response: {e}")
        print("Test 1 failed!")
        return False
    
    # Test 2: Generate a recipe using a structured prompt
    print("\nTest 2: Generate a recipe")
    recipe_prompt = """
    Create a recipe for a vegetarian pasta dish. 
    Format the response as a valid JSON object with these fields:
    - title: The name of the dish
    - ingredients: A list of ingredients needed
    - instructions: A list of step-by-step instructions
    - prep_time_minutes: Preparation time in minutes (integer)
    - cook_time_minutes: Cooking time in minutes (integer)
    - difficulty: Difficulty level (easy, medium, or hard)
    
    JSON:
    """
    
    # Create a proper LLMRequest for the recipe
    recipe_request = LLMRequest(
        user_prompt=recipe_prompt,
        llm_config=LLMModelConfig(
            provider=ProviderType.GROQ,
            model_name="llama3-70b-8192",  # Using Llama 3 70B model
            api_key=api_key
        ),
        max_tokens=1000
    )
    
    try:
        recipe_response = provider.generate(recipe_request)
        print("Raw response:")
        print(recipe_response.content)
        
        # Try to parse the JSON from the response
        try:
            # Find JSON in the response
            json_str = recipe_response.content
            if "```json" in recipe_response.content:
                json_str = recipe_response.content.split("```json")[1].split("```")[0].strip()
            elif "```" in recipe_response.content:
                json_str = recipe_response.content.split("```")[1].split("```")[0].strip()
            
            recipe_data = json.loads(json_str)
            recipe = Recipe(**recipe_data)
            print("\nParsed Recipe:")
            print(f"Title: {recipe.title}")
            print(f"Ingredients: {len(recipe.ingredients)} items")
            print(f"Instructions: {len(recipe.instructions)} steps")
            print(f"Prep time: {recipe.prep_time_minutes} minutes")
            print(f"Cook time: {recipe.cook_time_minutes} minutes")
            print(f"Difficulty: {recipe.difficulty}")
            print("Test 2 passed!")
        except Exception as e:
            print(f"Error parsing JSON: {e}")
            print("Test 2 failed!")
            return False
    except Exception as e:
        print(f"Error generating response: {e}")
        print("Test 2 failed!")
        return False
    
    print("\nAll tests passed!")
    return True

if __name__ == "__main__":
    test_groq_structured_simple() 