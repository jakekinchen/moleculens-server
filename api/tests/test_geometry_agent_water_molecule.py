"""
Test the geometry agent with different models using a water molecule prompt.
"""

import os
import asyncio
import json
import sys
from dotenv import load_dotenv
from api.agent_management.agents.geometry_agent import GeometryAgent
from api.agent_management.llm_service import LLMService, LLMModelConfig, LLMRequest, ProviderType

# Load environment variables
load_dotenv()

# Ensure we have API keys
os.environ["OPENAI_API_KEY"] = os.environ.get("OPENAI_API_KEY", "fake-key-for-testing")
os.environ["GROQ_API_KEY"] = os.environ.get("GROQ_API_KEY", "fake-key-for-testing")

# Define the test prompt for a water molecule
test_prompt = "Visualize a water molecule (H2O) with appropriate atomic sizes and bond angles."

# List of models to test with provider information (skipping Anthropic as provider is missing)
models_to_test = [
    {"name": "gpt-4o", "provider": ProviderType.OPENAI},
    {"name": "llama3-70b-8192", "provider": ProviderType.GROQ},
    {"name": "deepseek-r1-distill-llama-70b", "provider": ProviderType.GROQ}
]

async def test_model(model_info):
    """Test a specific model with the water molecule prompt."""
    model_name = model_info["name"]
    provider = model_info["provider"]
    
    print(f"\n\n{'=' * 50}")
    print(f"Testing model: {model_name} with provider: {provider}")
    print(f"{'=' * 50}\n")
    
    try:
        # Create LLM service with the model
        llm_service = LLMService(
            config=LLMModelConfig(
                provider=provider,
                model_name=model_name,
                api_key=os.getenv(f"{provider.upper()}_API_KEY")
            )
        )
        
        # Create a geometry agent with the service
        geometry_agent = GeometryAgent(llm_service=llm_service)
        
        # Generate geometry for water molecule
        result = geometry_agent.get_geometry_snippet(test_prompt)
        
        # Print the result
        print(f"Result from {model_name}:\n")
        print(result)
        
        # Return the result
        return {
            "model": model_name,
            "provider": provider,
            "success": True,
            "result": result
        }
    except Exception as e:
        print(f"Error testing {model_name}: {str(e)}")
        return {
            "model": model_name,
            "provider": provider,
            "success": False,
            "error": str(e)
        }

async def main():
    """Test all models and generate a report."""
    results = []
    
    for model_info in models_to_test:
        try:
            result = await test_model(model_info)
            results.append(result)
        except Exception as e:
            print(f"Failed to test {model_info['name']}: {str(e)}")
            results.append({
                "model": model_info["name"],
                "provider": model_info["provider"],
                "success": False,
                "error": str(e)
            })
    
    # Generate report
    report = {
        "test_prompt": test_prompt,
        "results": results
    }
    
    # Save report to file
    with open("geometry_agent_water_molecule_test_results.json", "w") as f:
        json.dump(report, f, indent=2)
    
    # Generate a readable report
    with open("geometry_agent_water_molecule_test_report.md", "w") as f:
        f.write("# Geometry Agent Test Results: Water Molecule\n\n")
        f.write(f"**Test prompt:** {test_prompt}\n\n")
        
        for result in results:
            model = result["model"]
            provider = result["provider"]
            f.write(f"## Model: {model} (Provider: {provider})\n\n")
            
            if result["success"]:
                f.write("**Status:** ✅ Success\n\n")
                f.write("**Generated code:**\n\n```javascript\n")
                f.write(result["result"])
                f.write("\n```\n\n")
            else:
                f.write("**Status:** ❌ Failed\n\n")
                f.write(f"**Error:** {result.get('error', 'Unknown error')}\n\n")
        
        f.write("\n## Summary\n\n")
        successful = sum(1 for r in results if r["success"])
        f.write(f"Successful models: {successful}/{len(models_to_test)}\n")
        
    print("\n\nTest completed. Results saved to:")
    print("- geometry_agent_water_molecule_test_results.json")
    print("- geometry_agent_water_molecule_test_report.md")

if __name__ == "__main__":
    asyncio.run(main())