"""
Test script to verify the Groq provider's JSON parsing capabilities
"""

import os
import sys
import json
from typing import List, Dict, Optional
from pydantic import BaseModel, Field
from dotenv import load_dotenv

# Add the parent directory to the path so we can import the modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api.agent_management.providers.groq_provider import GroqProvider
from api.agent_management.llm_service import LLMRequest, LLMModelConfig, ProviderType, StructuredLLMRequest

# Load environment variables
load_dotenv()

# Define test models for structured output
class SimpleModel(BaseModel):
    name: str
    value: int
    description: str

class NestedModel(BaseModel):
    title: str
    items: List[SimpleModel]
    metadata: Dict[str, str]

class ComplexModel(BaseModel):
    id: str
    nested: NestedModel
    tags: List[str]
    active: bool
    score: float = Field(ge=0.0, le=1.0)

def test_extract_json_from_text():
    """Test the _extract_json_from_text method directly"""
    print("\n" + "="*50)
    print("Testing _extract_json_from_text Method")
    print("="*50)
    
    # Create a Groq provider instance
    provider = GroqProvider()
    
    # Test cases with different JSON formats
    test_cases = [
        # Case 1: JSON in code block with json tag
        {
            "input": """Here's the JSON you requested:
```json
{
  "name": "Test Item",
  "value": 42,
  "description": "This is a test"
}
```
I hope this helps!""",
            "expected": '{\n  "name": "Test Item",\n  "value": 42,\n  "description": "This is a test"\n}'
        },
        
        # Case 2: JSON in code block without json tag
        {
            "input": """Here's the JSON:
```
{
  "name": "Test Item",
  "value": 42,
  "description": "This is a test"
}
```""",
            "expected": '{\n  "name": "Test Item",\n  "value": 42,\n  "description": "This is a test"\n}'
        },
        
        # Case 3: JSON directly in text
        {
            "input": """{"name": "Test Item", "value": 42, "description": "This is a test"}""",
            "expected": '{"name": "Test Item", "value": 42, "description": "This is a test"}'
        },
        
        # Case 4: JSON with prefix text
        {
            "input": """Here is a JSON object that conforms to the schema:
{
  "name": "Test Item",
  "value": 42,
  "description": "This is a test"
}""",
            "expected": '{\n  "name": "Test Item",\n  "value": 42,\n  "description": "This is a test"\n}'
        }
    ]
    
    # Run the tests
    for i, test_case in enumerate(test_cases):
        print(f"\nTest Case {i+1}:")
        result = provider._extract_json_from_text(test_case["input"])
        
        # Compare the result with the expected output
        try:
            # Parse both as JSON to compare the structure rather than exact string
            result_json = json.loads(result)
            expected_json = json.loads(test_case["expected"])
            
            if result_json == expected_json:
                print(f"✅ JSON extracted correctly")
            else:
                print(f"❌ JSON extraction failed")
                print(f"Expected: {test_case['expected']}")
                print(f"Got: {result}")
                return False
        except json.JSONDecodeError as e:
            print(f"❌ JSON extraction failed: {str(e)}")
            print(f"Expected: {test_case['expected']}")
            print(f"Got: {result}")
            return False
    
    return True

def test_fix_json_content():
    """Test the _fix_json_content method directly"""
    print("\n" + "="*50)
    print("Testing _fix_json_content Method")
    print("="*50)
    
    # Create a Groq provider instance
    provider = GroqProvider()
    
    # Test cases with different JSON issues
    test_cases = [
        # Case 1: JSON with Math.PI
        {
            "input": '{"angle": Math.PI, "half_angle": Math.PI / 2}',
            "expected": '{"angle": 3.141592653589793, "half_angle": 1.5707963267948966}'
        },
        
        # Case 2: JSON with single quotes
        {
            "input": "{'name': 'Test Item', 'value': 42}",
            "expected": '{"name": "Test Item", "value": 42}'
        },
        
        # Case 3: JSON with trailing commas
        {
            "input": '{"items": ["a", "b", "c",], "values": {"x": 1, "y": 2,}}',
            "expected": '{"items": ["a", "b", "c"], "values": {"x": 1, "y": 2}}'
        },
        
        # Case 4: JSON with Math functions
        {
            "input": '{"value": Math.sin(0), "random": Math.random()}',
            "expected": '{"value": 0, "random": 0}'
        }
    ]
    
    # Run the tests
    for i, test_case in enumerate(test_cases):
        print(f"\nTest Case {i+1}:")
        result = provider._fix_json_content(test_case["input"])
        
        # Compare the result with the expected output
        try:
            # Parse both as JSON to compare the structure rather than exact string
            result_json = json.loads(result)
            expected_json = json.loads(test_case["expected"])
            
            if result_json == expected_json:
                print(f"✅ JSON fixed correctly")
            else:
                print(f"❌ JSON fixing failed")
                print(f"Expected: {test_case['expected']}")
                print(f"Got: {result}")
                return False
        except json.JSONDecodeError as e:
            print(f"❌ JSON fixing failed: {str(e)}")
            print(f"Expected: {test_case['expected']}")
            print(f"Got: {result}")
            return False
    
    return True

def test_structured_output():
    """Test structured output generation with the Groq provider"""
    print("\n" + "="*50)
    print("Testing Structured Output Generation")
    print("="*50)
    
    # Create a Groq provider instance
    provider = GroqProvider()
    
    # Create a model config
    model_config = LLMModelConfig(
        provider=ProviderType.GROQ,
        model_name="llama3-8b-8192"
    )
    
    # Test with SimpleModel
    print("\nTest with SimpleModel:")
    simple_request = StructuredLLMRequest(
        user_prompt="Create a test item with name 'Test Item', value 42, and a short description.",
        system_prompt="You are a helpful assistant that generates structured data.",
        response_model=SimpleModel,
        llm_config=model_config
    )
    
    try:
        simple_result = provider.generate_structured(simple_request)
        print(f"✅ Successfully generated SimpleModel")
        print(f"Result: {simple_result}")
    except Exception as e:
        print(f"❌ Failed to generate SimpleModel: {str(e)}")
        return False
    
    # Test with NestedModel
    print("\nTest with NestedModel:")
    nested_request = StructuredLLMRequest(
        user_prompt="""Create a list with title 'Test List', 2 items, and metadata with 'created_by' and 'date' fields.

IMPORTANT: Your response must be a valid JSON object with the following structure:
{
  "title": "Test List",
  "items": [
    {
      "name": "Item 1",
      "value": 1,
      "description": "Description 1"
    },
    {
      "name": "Item 2",
      "value": 2,
      "description": "Description 2"
    }
  ],
  "metadata": {
    "created_by": "Author Name",
    "date": "2023-01-01"
  }
}

Do NOT include the schema definition in your response. Only return a valid instance of the model.""",
        system_prompt="You are a helpful assistant that generates structured data. Always return valid JSON that matches the requested structure exactly.",
        response_model=NestedModel,
        llm_config=model_config
    )
    
    try:
        nested_result = provider.generate_structured(nested_request)
        print(f"✅ Successfully generated NestedModel")
        print(f"Result title: {nested_result.title}")
        print(f"Result items count: {len(nested_result.items)}")
        print(f"Result metadata: {nested_result.metadata}")
    except Exception as e:
        print(f"❌ Failed to generate NestedModel: {str(e)}")
        return False
    
    # Test with ComplexModel
    print("\nTest with ComplexModel:")
    complex_request = StructuredLLMRequest(
        user_prompt="""Create a complex object with an ID, a nested structure, tags, active status, and a score between 0 and 1.

IMPORTANT: Your response must be a valid JSON object with the following structure:
{
  "id": "abc123",
  "nested": {
    "title": "Nested Title",
    "items": [
      {
        "name": "Item 1",
        "value": 1,
        "description": "Description 1"
      },
      {
        "name": "Item 2",
        "value": 2,
        "description": "Description 2"
      }
    ],
    "metadata": {
      "created_by": "Author Name",
      "date": "2023-01-01"
    }
  },
  "tags": ["tag1", "tag2", "tag3"],
  "active": true,
  "score": 0.75
}

Do NOT include the schema definition in your response. Only return a valid instance of the model.""",
        system_prompt="You are a helpful assistant that generates structured data. Always return valid JSON that matches the requested structure exactly.",
        response_model=ComplexModel,
        llm_config=model_config
    )
    
    try:
        complex_result = provider.generate_structured(complex_request)
        print(f"✅ Successfully generated ComplexModel")
        print(f"Result ID: {complex_result.id}")
        print(f"Result nested title: {complex_result.nested.title}")
        print(f"Result tags count: {len(complex_result.tags)}")
        print(f"Result active: {complex_result.active}")
        print(f"Result score: {complex_result.score}")
    except Exception as e:
        print(f"❌ Failed to generate ComplexModel: {str(e)}")
        return False
    
    return True

def test_groq_json_parsing():
    """Test the Groq provider's JSON parsing capabilities"""
    print("\n" + "="*50)
    print("Testing Groq Provider JSON Parsing")
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
    
    # Run the tests
    test_results = {}
    
    # Test _extract_json_from_text
    try:
        extract_result = test_extract_json_from_text()
        test_results["_extract_json_from_text"] = "✅ Passed" if extract_result else "❌ Failed"
    except Exception as e:
        print(f"Error testing _extract_json_from_text: {str(e)}")
        test_results["_extract_json_from_text"] = f"❌ Failed: {str(e)}"
    
    # Test _fix_json_content
    try:
        fix_result = test_fix_json_content()
        test_results["_fix_json_content"] = "✅ Passed" if fix_result else "❌ Failed"
    except Exception as e:
        print(f"Error testing _fix_json_content: {str(e)}")
        test_results["_fix_json_content"] = f"❌ Failed: {str(e)}"
    
    # Test structured output generation
    try:
        structured_result = test_structured_output()
        test_results["structured_output"] = "✅ Passed" if structured_result else "❌ Failed"
    except Exception as e:
        print(f"Error testing structured_output: {str(e)}")
        test_results["structured_output"] = f"❌ Failed: {str(e)}"
    
    # Print summary of results
    print("\n" + "="*50)
    print("SUMMARY OF JSON PARSING TEST RESULTS")
    print("="*50)
    for test, result in test_results.items():
        print(f"{test}: {result}")
    
    # Overall result
    if all("Passed" in result for result in test_results.values()):
        print("\n✅ All JSON parsing tests passed with Groq provider")
        return True
    else:
        print("\n❌ Some JSON parsing tests failed with Groq provider")
        return False

if __name__ == "__main__":
    success = test_groq_json_parsing()
    if success:
        print("\n✅ All tests passed - Groq provider JSON parsing is working correctly")
    else:
        print("\n❌ Tests failed - Groq provider JSON parsing has issues") 