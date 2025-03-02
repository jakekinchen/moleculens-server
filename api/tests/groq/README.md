# Groq Provider Tests

This directory contains tests for the Groq provider and its models.

## Test Files

### Core Tests
- `test_groq_model_registry.py`: Tests that all models are correctly registered in the model registry
- `test_groq_new_models.py`: Tests basic functionality of all models

### Capability Tests
- `test_groq_code_model.py`: Tests code generation capabilities of the DeepSeek R1 Distill Llama 70B model
- `test_groq_reasoning_model.py`: Tests reasoning capabilities of the Llama 3.3 70B SpecDec model
- `test_groq_structured_output.py`: Tests structured output generation with the Qwen 2.5 32B model
- `test_groq_vision_model.py`: Tests vision capabilities of the Llama 3.2 90B Vision model

### Provider Tests
- `test_groq_simple.py`: Tests basic provider functionality
- `test_groq_service.py`: Tests provider through LLMService
- `test_groq_structured.py`: Tests structured output generation
- `test_groq_all_models.py`: Tests all response models
- `test_groq_structured_simple.py`: Tests basic structured output without instructor
- `test_groq_api_key_management.py`: Tests API key management
- `test_groq_provider_implementation.py`: Tests provider implementation
- `test_groq_provider_registration.py`: Tests provider registration
- `test_groq_json_parsing.py`: Tests JSON parsing
- `test_groq_instructor.py`: Tests instructor integration

### Integration Tests
- `test_groq_agents.py`: Tests agent compatibility

### Test Runner
- `run_all_tests.py`: Runs all tests and reports results

## Running Tests

### Running Individual Tests

To run an individual test, use:

```bash
python api/tests/groq/test_groq_model_registry.py
```

### Running All Tests

To run all tests, use:

```bash
python api/tests/groq/run_all_tests.py
```

### Running Tests by Category

To run tests of a specific category, use:

```bash
python api/tests/groq/run_all_tests.py --category core
```

Available categories:
- `core`: Core functionality tests
- `capability`: Model capability tests
- `provider`: Provider implementation tests
- `integration`: Integration tests

### Generating a Test Report

To run all tests and generate a comprehensive report, use:

```bash
python api/tests/groq/run_tests_and_report.py
```

This will:
1. Run all tests in each category
2. Generate a Markdown report with test results
3. Save the report in the `reports` directory
4. Attempt to open the report automatically

## Test Requirements

- The Groq API key must be set in the environment variable `GROQ_API_KEY`
- For the vision test, an image file is required in the `test_images` directory

## Test Results

The test results are documented in:

- `api/agent_management/providers/GROQ_TESTING_SUMMARY.md`: Summary of test results
- `api/agent_management/providers/GROQ_TESTING_FINAL_REPORT.md`: Comprehensive test report

## Known Issues

- The vision test may not work correctly as the Groq API might not properly support image input yet
- Some tests may be truncated due to token limitations 