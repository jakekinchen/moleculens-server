#!/usr/bin/env python3
"""
Main test runner script for Groq tests.

This script runs all the Groq tests in sequence and reports the results.
"""

import os
import sys
import importlib
from typing import List, Dict, Any, Tuple, Optional
import argparse

# Add the parent directory to the path so we can import the modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

# Import test modules
try:
    import api.tests.groq.test_groq_model_registry as test_groq_model_registry
    import api.tests.groq.test_groq_new_models as test_groq_new_models
    import api.tests.groq.test_groq_code_model as test_groq_code_model
    import api.tests.groq.test_groq_reasoning_model as test_groq_reasoning_model
    import api.tests.groq.test_groq_structured_output as test_groq_structured_output
    import api.tests.groq.test_groq_vision_model as test_groq_vision_model

    # Additional test modules
    import api.tests.groq.test_groq_simple as test_groq_simple
    import api.tests.groq.test_groq_service as test_groq_service
    import api.tests.groq.test_groq_structured as test_groq_structured
    import api.tests.groq.test_groq_all_models as test_groq_all_models
    import api.tests.groq.test_groq_structured_simple as test_groq_structured_simple
    import api.tests.groq.test_groq_api_key_management as test_groq_api_key_management
    import api.tests.groq.test_groq_provider_implementation as test_groq_provider_implementation
    import api.tests.groq.test_groq_provider_registration as test_groq_provider_registration
    import api.tests.groq.test_groq_json_parsing as test_groq_json_parsing
    import api.tests.groq.test_groq_instructor as test_groq_instructor
    import api.tests.groq.test_groq_agents as test_groq_agents
except ImportError as e:
    print(f"Warning: Some test modules could not be imported: {e}")

def run_test(test_module, test_function_name: str) -> Tuple[bool, str]:
    """Run a test function from a module and return the result."""
    try:
        test_function = getattr(test_module, test_function_name)
        result = test_function()
        return result, ""
    except Exception as e:
        return False, str(e)

def run_all_tests(category: Optional[str] = None) -> Dict[str, Dict[str, Any]]:
    """Run all Groq tests and return the results."""
    tests = [
        # Core tests
        {
            "name": "Model Registry Test",
            "module": test_groq_model_registry,
            "function": "test_groq_model_registry",
            "category": "core"
        },
        {
            "name": "New Models Test",
            "module": test_groq_new_models,
            "function": "test_all_new_models",
            "category": "core"
        },
        {
            "name": "Code Generation Test",
            "module": test_groq_code_model,
            "function": "test_code_generation",
            "category": "capability"
        },
        {
            "name": "Reasoning Test",
            "module": test_groq_reasoning_model,
            "function": "test_reasoning_capabilities",
            "category": "capability"
        },
        {
            "name": "Structured Output Test",
            "module": test_groq_structured_output,
            "function": "test_structured_output",
            "category": "capability"
        },
        {
            "name": "Vision Test",
            "module": test_groq_vision_model,
            "function": "test_vision_model",
            "category": "capability"
        }
    ]
    
    # Add additional tests if the modules are available
    additional_tests = []
    
    if 'test_groq_simple' in globals():
        additional_tests.append({
            "name": "Simple Provider Test",
            "module": test_groq_simple,
            "function": "test_groq_simple",
            "category": "provider"
        })
    
    if 'test_groq_service' in globals():
        additional_tests.append({
            "name": "Service Integration Test",
            "module": test_groq_service,
            "function": "test_groq_service",
            "category": "provider"
        })
    
    if 'test_groq_structured' in globals():
        additional_tests.append({
            "name": "Structured Output Provider Test",
            "module": test_groq_structured,
            "function": "test_groq_structured",
            "category": "provider"
        })
    
    if 'test_groq_all_models' in globals():
        additional_tests.append({
            "name": "All Models Test",
            "module": test_groq_all_models,
            "function": "test_groq_all_models",
            "category": "provider"
        })
    
    if 'test_groq_structured_simple' in globals():
        additional_tests.append({
            "name": "Simple Structured Output Test",
            "module": test_groq_structured_simple,
            "function": "test_groq_structured_simple",
            "category": "provider"
        })
    
    if 'test_groq_api_key_management' in globals():
        additional_tests.append({
            "name": "API Key Management Test",
            "module": test_groq_api_key_management,
            "function": "test_groq_api_key_management",
            "category": "provider"
        })
    
    if 'test_groq_provider_implementation' in globals():
        additional_tests.append({
            "name": "Provider Implementation Test",
            "module": test_groq_provider_implementation,
            "function": "test_groq_provider_implementation",
            "category": "provider"
        })
    
    if 'test_groq_provider_registration' in globals():
        additional_tests.append({
            "name": "Provider Registration Test",
            "module": test_groq_provider_registration,
            "function": "test_groq_provider_registration",
            "category": "provider"
        })
    
    if 'test_groq_json_parsing' in globals():
        additional_tests.append({
            "name": "JSON Parsing Test",
            "module": test_groq_json_parsing,
            "function": "test_groq_json_parsing",
            "category": "provider"
        })
    
    if 'test_groq_instructor' in globals():
        additional_tests.append({
            "name": "Instructor Test",
            "module": test_groq_instructor,
            "function": "test_groq_instructor",
            "category": "provider"
        })
    
    if 'test_groq_agents' in globals():
        additional_tests.append({
            "name": "Agents Test",
            "module": test_groq_agents,
            "function": "test_groq_agents",
            "category": "integration"
        })
    
    # Add additional tests to the main tests list
    tests.extend(additional_tests)
    
    # Filter tests by category if specified
    if category:
        tests = [test for test in tests if test["category"] == category]
    
    results = {}
    
    print("\n=== Running Groq Tests ===\n")
    if category:
        print(f"Category: {category}\n")
    
    for test in tests:
        print(f"Running {test['name']}...")
        success, error = run_test(test["module"], test["function"])
        results[test["name"]] = {
            "success": success,
            "error": error,
            "category": test["category"]
        }
        if success:
            print(f"✅ {test['name']} passed!\n")
        else:
            print(f"❌ {test['name']} failed: {error}\n")
    
    # Print summary
    print("\n=== Test Summary ===\n")
    passed = sum(1 for test in results.values() if test["success"])
    total = len(results)
    if total > 0:
        print(f"Passed: {passed}/{total} tests ({passed/total*100:.1f}%)\n")
    else:
        print("No tests were run.\n")
    
    # Group results by category
    categories = {}
    for name, result in results.items():
        category = result["category"]
        if category not in categories:
            categories[category] = {"passed": 0, "total": 0}
        categories[category]["total"] += 1
        if result["success"]:
            categories[category]["passed"] += 1
    
    # Print results by category
    for category, stats in categories.items():
        print(f"{category.capitalize()} Tests: {stats['passed']}/{stats['total']} passed")
    
    print("\nDetailed Results:")
    for name, result in results.items():
        status = "✅ Passed" if result["success"] else f"❌ Failed: {result['error']}"
        print(f"{name} ({result['category']}): {status}")
    
    return results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Groq tests")
    parser.add_argument("--category", type=str, help="Run tests of a specific category (core, capability, provider, integration)")
    args = parser.parse_args()
    
    run_all_tests(args.category) 