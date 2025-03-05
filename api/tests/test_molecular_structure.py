#!/usr/bin/env python3
"""
Test script to evaluate the molecular_structure method with different LLM models.
"""

import os
import sys
import json
import time
from typing import Dict, Any, List, Optional
from datetime import datetime

# Add the parent directory to the path so we can import the modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api.agent_management.model_config import ModelRegistry
from api.agent_management.llm_service import LLMService, LLMModelConfig, ProviderType
from api.agent_management.agents.domain_bool_agent import DomainValidator
from api.agent_management.providers.deepseek_utils import is_deepseek_model

# Test prompts
TEST_PROMPTS = [
    "Draw a molecule of water",
    "Show the structure of aspirin",
    "Create a methane molecule"
]

def is_valid_smiles(smiles: str) -> bool:
    """
    Very basic check if a string looks like a SMILES string.
    This is not a comprehensive validation, just a basic check.
    """
    # Check if the string contains only valid SMILES characters
    valid_chars = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789()[]{}@+-=\\/~#$%^&*.")
    return all(c in valid_chars for c in smiles) and len(smiles) > 0

def test_model(model_name: str, prompt: str) -> Dict[str, Any]:
    """
    Test a model with a prompt and return the results.
    """
    result = {
        "model_name": model_name,
        "prompt": prompt,
        "success": False,
        "is_deepseek": False,
        "response": None,
        "error": None,
        "is_valid_smiles": False,
        "elapsed_time": 0
    }
    
    try:
        # Get model info
        model_info = ModelRegistry.create_instance(model_name)
        result["provider"] = model_info.provider.value
        result["is_deepseek"] = is_deepseek_model(model_info.name)
        
        # Create LLM service
        llm_service = LLMService(
            config=LLMModelConfig(
                provider=model_info.provider,
                model_name=model_info.name,
                api_key=None  # API key will be loaded from environment variables
            )
        )
        
        # Create domain validator
        domain_validator = DomainValidator(llm_service)
        
        # Time the request
        start_time = time.time()
        
        # Get molecular structure
        response = domain_validator.molecular_structure(prompt)
        
        # Calculate elapsed time
        result["elapsed_time"] = time.time() - start_time
        
        # Check if response is valid
        if hasattr(response, "molecular_structure"):
            result["success"] = True
            result["response"] = response.molecular_structure
            result["is_valid_smiles"] = is_valid_smiles(response.molecular_structure)
        else:
            result["success"] = False
            result["response"] = str(response)
            result["error"] = "Response does not have molecular_structure attribute"
            
    except Exception as e:
        result["success"] = False
        result["error"] = str(e)
    
    return result

def run_tests() -> List[Dict[str, Any]]:
    """
    Run tests for all models and prompts.
    """
    results = []
    
    # Get all models
    models = ModelRegistry.list_models()
    
    # Check environment variables for API keys
    missing_keys = []
    if not os.getenv("OPENAI_API_KEY"):
        missing_keys.append("OPENAI_API_KEY")
    if not os.getenv("ANTHROPIC_API_KEY"):
        missing_keys.append("ANTHROPIC_API_KEY")
    if not os.getenv("GROQ_API_KEY"):
        missing_keys.append("GROQ_API_KEY")
    
    if missing_keys:
        print(f"Warning: The following API keys are missing: {', '.join(missing_keys)}")
        print("Some tests may fail due to missing API keys.")
    
    # Test each model with each prompt
    for model_name in models:
        model_info = ModelRegistry.create_instance(model_name)
        print(f"Testing model: {model_name} ({model_info.display_name})")
        
        for prompt in TEST_PROMPTS:
            print(f"  Prompt: {prompt}")
            result = test_model(model_name, prompt)
            results.append(result)
            
            # Print result
            if result["success"]:
                if result["is_valid_smiles"]:
                    print(f"    ✅ Success: {result['response']}")
                else:
                    print(f"    ⚠️ Invalid SMILES: {result['response']}")
            else:
                print(f"    ❌ Error: {result['error']}")
            
            # Sleep to avoid rate limiting
            time.sleep(1)
    
    return results

def generate_report(results: List[Dict[str, Any]]) -> str:
    """
    Generate a report from the test results.
    """
    # Group results by model
    model_results = {}
    for result in results:
        model_name = result["model_name"]
        if model_name not in model_results:
            model_results[model_name] = []
        model_results[model_name].append(result)
    
    # Generate report
    report = "# Molecular Structure Test Report\n\n"
    report += f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"
    
    # Summary
    total_tests = len(results)
    successful_tests = sum(1 for r in results if r["success"])
    valid_smiles_tests = sum(1 for r in results if r["is_valid_smiles"])
    
    report += f"## Summary\n\n"
    report += f"- Total tests: {total_tests}\n"
    report += f"- Successful tests: {successful_tests} ({successful_tests/total_tests*100:.1f}%)\n"
    report += f"- Valid SMILES: {valid_smiles_tests} ({valid_smiles_tests/total_tests*100:.1f}%)\n\n"
    
    # Results by model
    report += f"## Results by Model\n\n"
    
    for model_name, model_results_list in model_results.items():
        model_info = ModelRegistry.create_instance(model_name)
        
        # Calculate statistics
        total_model_tests = len(model_results_list)
        successful_model_tests = sum(1 for r in model_results_list if r["success"])
        valid_smiles_model_tests = sum(1 for r in model_results_list if r["is_valid_smiles"])
        
        success_rate = successful_model_tests / total_model_tests * 100 if total_model_tests > 0 else 0
        valid_smiles_rate = valid_smiles_model_tests / total_model_tests * 100 if total_model_tests > 0 else 0
        
        # Add model section
        report += f"### {model_info.display_name} ({model_name})\n\n"
        report += f"- Provider: {model_results_list[0]['provider']}\n"
        report += f"- Is DeepSeek model: {'Yes' if model_results_list[0]['is_deepseek'] else 'No'}\n"
        report += f"- Success rate: {success_rate:.1f}%\n"
        report += f"- Valid SMILES rate: {valid_smiles_rate:.1f}%\n\n"
        
        # Add results table
        report += "| Prompt | Success | Response | Valid SMILES | Time (s) |\n"
        report += "|--------|---------|----------|--------------|----------|\n"
        
        for result in model_results_list:
            success = "✅" if result["success"] else "❌"
            valid_smiles = "✅" if result["is_valid_smiles"] else "❌"
            response = result["response"] if result["response"] else result["error"]
            response = response.replace("\n", " ")[:50]  # Truncate long responses
            
            report += f"| {result['prompt']} | {success} | {response} | {valid_smiles} | {result['elapsed_time']:.2f} |\n"
        
        report += "\n"
    
    # Conclusion
    report += "## Conclusion\n\n"
    
    # Find best model
    best_model = None
    best_valid_rate = -1
    
    for model_name, model_results_list in model_results.items():
        valid_smiles_model_tests = sum(1 for r in model_results_list if r["is_valid_smiles"])
        valid_smiles_rate = valid_smiles_model_tests / len(model_results_list) * 100
        
        if valid_smiles_rate > best_valid_rate:
            best_valid_rate = valid_smiles_rate
            best_model = model_name
    
    if best_model:
        model_info = ModelRegistry.create_instance(best_model)
        report += f"The best performing model was **{model_info.display_name}** with a valid SMILES rate of {best_valid_rate:.1f}%.\n\n"
    
    # DeepSeek models analysis
    deepseek_results = [r for r in results if r["is_deepseek"]]
    if deepseek_results:
        deepseek_success = sum(1 for r in deepseek_results if r["success"])
        deepseek_valid = sum(1 for r in deepseek_results if r["is_valid_smiles"])
        
        report += f"### DeepSeek Models Analysis\n\n"
        report += f"- Total DeepSeek tests: {len(deepseek_results)}\n"
        report += f"- Successful tests: {deepseek_success} ({deepseek_success/len(deepseek_results)*100:.1f}%)\n"
        report += f"- Valid SMILES: {deepseek_valid} ({deepseek_valid/len(deepseek_results)*100:.1f}%)\n\n"
        
        report += "This indicates that the DeepSeek extraction utilities are "
        if deepseek_valid / len(deepseek_results) > 0.5:
            report += "working effectively for structured responses.\n"
        else:
            report += "not working effectively for structured responses and may need improvement.\n"
    
    return report

if __name__ == "__main__":
    print("Running molecular structure tests...")
    results = run_tests()
    
    # Generate report
    report = generate_report(results)
    
    # Save report to file
    report_file = "molecular_structure_test_report.md"
    with open(report_file, "w") as f:
        f.write(report)
    
    print(f"\nReport saved to {report_file}")
    
    # Also save raw results as JSON for further analysis
    results_file = "molecular_structure_test_results.json"
    with open(results_file, "w") as f:
        # Convert results to serializable format
        serializable_results = []
        for result in results:
            serializable_result = result.copy()
            if hasattr(serializable_result["response"], "model_dump"):
                serializable_result["response"] = serializable_result["response"].model_dump()
            serializable_results.append(serializable_result)
        
        json.dump(serializable_results, f, indent=2)
    
    print(f"Raw results saved to {results_file}") 