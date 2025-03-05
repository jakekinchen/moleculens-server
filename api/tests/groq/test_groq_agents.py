"""
Test script to verify the Groq provider works with agents
"""

import os
import sys
from typing import Dict, Any
from dotenv import load_dotenv

# Add the parent directory to the path so we can import the modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api.agent_management.llm_service import LLMService, LLMModelConfig, ProviderType
from api.agent_management.agents.domain_bool_agent import DomainValidator
from api.agent_management.agents.geometry_agent import GeometryAgent

# Load environment variables
load_dotenv()

def test_domain_validator():
    """Test the DomainValidator agent with the Groq provider"""
    print("\n" + "="*50)
    print("Testing DomainValidator with Groq Provider")
    print("="*50)
    
    # Create an LLM service
    llm_service = LLMService()
    
    # Create a model config for Groq
    model_config = LLMModelConfig(
        provider=ProviderType.GROQ,
        model_name="llama3-8b-8192"
    )
    
    # Create a DomainValidator
    validator = DomainValidator(llm_service=llm_service)
    
    # Test with a scientific prompt
    scientific_prompt = "What is the molecular structure of ethanol?"
    print(f"\nTesting with scientific prompt: '{scientific_prompt}'")
    
    try:
        result = validator.is_molecular(scientific_prompt)
        print(f"Is scientific: {result.is_true}")
        print(f"Confidence: {result.confidence}")
        print(f"Reasoning: {result.reasoning[:100]}...")
        print("✅ DomainValidator test passed")
    except Exception as e:
        print(f"❌ DomainValidator test failed: {str(e)}")
        return False
    
    # Test with a non-scientific prompt
    non_scientific_prompt = "What is the best recipe for chocolate chip cookies?"
    print(f"\nTesting with non-scientific prompt: '{non_scientific_prompt}'")
    
    try:
        result = validator.is_molecular(non_scientific_prompt)
        print(f"Is scientific: {result.is_true}")
        print(f"Confidence: {result.confidence}")
        print(f"Reasoning: {result.reasoning[:100]}...")
        print("✅ DomainValidator test passed")
    except Exception as e:
        print(f"❌ DomainValidator test failed: {str(e)}")
        return False
    
    return True

def test_geometry_agent():
    """Test the GeometryAgent with the Groq provider"""
    print("\n" + "="*50)
    print("Testing GeometryAgent with Groq Provider")
    print("="*50)
    
    # Create an LLM service
    llm_service = LLMService()
    
    # Create a model config for Groq
    model_config = LLMModelConfig(
        provider=ProviderType.GROQ,
        model_name="llama3-8b-8192"
    )
    
    # Create a GeometryAgent
    geometry_agent = GeometryAgent(llm_service=llm_service)
    
    # Test with a simple geometry prompt
    geometry_prompt = "Create a sphere with radius 5"
    print(f"\nTesting with geometry prompt: '{geometry_prompt}'")
    
    try:
        geometry_code = geometry_agent.get_geometry_snippet(geometry_prompt)
        print(f"Generated code length: {len(geometry_code)}")
        print("Sample of generated code:")
        print(geometry_code[:200] + "..." if len(geometry_code) > 200 else geometry_code)
        print("✅ GeometryAgent test passed")
    except Exception as e:
        print(f"❌ GeometryAgent test failed: {str(e)}")
        return False
    
    return True

def test_groq_agents():
    """Test the Groq provider with agents"""
    print("\n" + "="*50)
    print("Testing Groq Provider with Agents")
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
    
    # Test DomainValidator
    try:
        domain_result = test_domain_validator()
        test_results["DomainValidator"] = "✅ Passed" if domain_result else "❌ Failed"
    except Exception as e:
        print(f"Error testing DomainValidator: {str(e)}")
        test_results["DomainValidator"] = f"❌ Failed: {str(e)}"
    
    # Test GeometryAgent
    try:
        geometry_result = test_geometry_agent()
        test_results["GeometryAgent"] = "✅ Passed" if geometry_result else "❌ Failed"
    except Exception as e:
        print(f"Error testing GeometryAgent: {str(e)}")
        test_results["GeometryAgent"] = f"❌ Failed: {str(e)}"
    
    # Print summary of results
    print("\n" + "="*50)
    print("SUMMARY OF AGENT TEST RESULTS")
    print("="*50)
    for test, result in test_results.items():
        print(f"{test}: {result}")
    
    # Overall result
    if all("Passed" in result for result in test_results.values()):
        print("\n✅ All agent tests passed with Groq provider")
        return True
    else:
        print("\n❌ Some agent tests failed with Groq provider")
        return False

if __name__ == "__main__":
    success = test_groq_agents()
    if success:
        print("\n✅ All tests passed - Groq provider works with agents")
    else:
        print("\n❌ Tests failed - Groq provider has issues with agents") 