"""
Test script to demonstrate using the instructor library with Groq for structured outputs.
"""

import os
import sys
from typing import List, Dict, Optional
from pydantic import BaseModel, Field
from dotenv import load_dotenv

# Add the parent directory to the path so we can import the modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Load environment variables
load_dotenv()

# Define test models for structured output
class UserInfo(BaseModel):
    name: str
    age: int
    email: str

class Recipe(BaseModel):
    title: str
    ingredients: List[str]
    steps: List[str]
    prep_time_minutes: int
    cook_time_minutes: int

def test_groq_instructor():
    """Test structured output with Groq using the instructor library"""
    print("\n" + "="*50)
    print("Testing Groq with Instructor for Structured Output")
    print("="*50)
    
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
    
    # Check if instructor package is installed
    try:
        import instructor
        print("Instructor package is installed")
    except ImportError:
        print("Error: instructor package is not installed")
        print("Install with: pip install instructor")
        return False
    
    # Create a Groq client
    client = groq.Groq(api_key=os.getenv("GROQ_API_KEY"))
    
    # Patch the client with instructor
    instructor_client = instructor.from_groq(client, mode=instructor.Mode.JSON)
    
    # Test with UserInfo model
    print("\nTesting with UserInfo model:")
    user_text = """
    John Doe, a 35-year-old software engineer from New York, has been working with large language models for several years.
    His email address is johndoe@example.com.
    """
    
    try:
        user_info = instructor_client.chat.completions.create(
            model="llama3-8b-8192",
            response_model=UserInfo,
            messages=[
                {"role": "system", "content": "Your job is to extract user information from the given text."},
                {"role": "user", "content": user_text}
            ],
            temperature=0.7,
        )
        
        print(f"✅ Successfully extracted UserInfo")
        print(f"Name: {user_info.name}")
        print(f"Age: {user_info.age}")
        print(f"Email: {user_info.email}")
    except Exception as e:
        print(f"❌ Failed to extract UserInfo: {str(e)}")
        return False
    
    # Test with Recipe model
    print("\nTesting with Recipe model:")
    recipe_prompt = "Provide a recipe for chocolate chip cookies."
    
    try:
        recipe = instructor_client.chat.completions.create(
            model="llama3-8b-8192",
            response_model=Recipe,
            messages=[
                {"role": "system", "content": "You are a helpful cooking assistant that provides recipes."},
                {"role": "user", "content": recipe_prompt}
            ],
            temperature=0.7,
        )
        
        print(f"✅ Successfully generated Recipe")
        print(f"Title: {recipe.title}")
        print(f"Ingredients: {len(recipe.ingredients)} items")
        print(f"Steps: {len(recipe.steps)} steps")
        print(f"Prep time: {recipe.prep_time_minutes} minutes")
        print(f"Cook time: {recipe.cook_time_minutes} minutes")
    except Exception as e:
        print(f"❌ Failed to generate Recipe: {str(e)}")
        return False
    
    print("\n✅ All tests passed - Instructor with Groq is working correctly")
    return True

if __name__ == "__main__":
    test_groq_instructor() 