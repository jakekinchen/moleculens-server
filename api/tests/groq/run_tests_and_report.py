#!/usr/bin/env python3
"""
Script to run all Groq tests and generate a report.
"""

import os
import sys
import subprocess
import datetime
from pathlib import Path

# Add the parent directory to the path so we can import the modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

# Import the run_all_tests function
from api.tests.groq.run_all_tests import run_all_tests

def generate_report():
    """Run all tests and generate a report."""
    # Create a reports directory if it doesn't exist
    reports_dir = Path(__file__).parent / "reports"
    reports_dir.mkdir(exist_ok=True)
    
    # Generate a timestamp for the report filename
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    report_file = reports_dir / f"groq_test_report_{timestamp}.md"
    
    # Run all tests
    print("Running all Groq tests...")
    results = {}
    
    # Run tests by category
    categories = ["core", "capability", "provider", "integration"]
    for category in categories:
        print(f"\nRunning {category} tests...")
        category_results = run_all_tests(category)
        results[category] = category_results
    
    # Generate the report
    print(f"\nGenerating report: {report_file}")
    with open(report_file, "w") as f:
        f.write("# Groq Tests Report\n\n")
        f.write(f"Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # Summary
        f.write("## Summary\n\n")
        total_tests = 0
        total_passed = 0
        
        for category, category_results in results.items():
            category_total = len(category_results)
            category_passed = sum(1 for test in category_results.values() if test["success"])
            total_tests += category_total
            total_passed += category_passed
            
            f.write(f"- **{category.capitalize()} Tests**: {category_passed}/{category_total} passed ")
            if category_total > 0:
                f.write(f"({category_passed/category_total*100:.1f}%)\n")
            else:
                f.write("(0.0%)\n")
        
        if total_tests > 0:
            f.write(f"\n**Overall**: {total_passed}/{total_tests} tests passed ({total_passed/total_tests*100:.1f}%)\n\n")
        else:
            f.write("\n**Overall**: No tests were run.\n\n")
        
        # Detailed results by category
        for category, category_results in results.items():
            f.write(f"## {category.capitalize()} Tests\n\n")
            
            if not category_results:
                f.write("No tests were run in this category.\n\n")
                continue
            
            for name, result in category_results.items():
                status = "✅ Passed" if result["success"] else f"❌ Failed: {result['error']}"
                f.write(f"### {name}\n\n")
                f.write(f"Status: {status}\n\n")
        
        # Recommendations
        f.write("## Recommendations\n\n")
        
        # Check if any tests failed
        failed_tests = []
        for category, category_results in results.items():
            for name, result in category_results.items():
                if not result["success"]:
                    failed_tests.append((name, category, result["error"]))
        
        if failed_tests:
            f.write("### Failed Tests\n\n")
            f.write("The following tests failed and should be investigated:\n\n")
            for name, category, error in failed_tests:
                f.write(f"- **{name}** ({category}): {error}\n")
            f.write("\n")
        
        # Check for vision test
        vision_test_found = False
        vision_test_passed = False
        for category, category_results in results.items():
            for name, result in category_results.items():
                if "Vision" in name:
                    vision_test_found = True
                    vision_test_passed = result["success"]
        
        if vision_test_found:
            f.write("### Vision Capabilities\n\n")
            if vision_test_passed:
                f.write("The vision test passed, but it's important to note that the current implementation doesn't fully support image input yet. The test is designed to handle this limitation gracefully.\n\n")
                f.write("Recommendations:\n\n")
                f.write("1. Monitor Groq's documentation for updates on the vision model's capabilities\n")
                f.write("2. Consider implementing a fallback to a text-only model if vision features are not critical\n")
                f.write("3. Investigate alternative methods for passing images to the API\n\n")
            else:
                f.write("The vision test failed. This suggests that there might be issues with the vision capabilities of the model or with the implementation.\n\n")
                f.write("Recommendations:\n\n")
                f.write("1. Check the Groq API documentation for the correct way to pass images to the model\n")
                f.write("2. Verify that the model supports vision capabilities\n")
                f.write("3. Consider using a different model for vision tasks\n\n")
        
        # General recommendations
        f.write("### General Recommendations\n\n")
        f.write("1. **Context Length Utilization**: Implement strategies to leverage the large context windows of these models (128K tokens for most models)\n")
        f.write("2. **Performance Monitoring**: Implement monitoring for response times and token usage\n")
        f.write("3. **Documentation**: Keep documentation updated as models evolve\n")
        f.write("4. **API Key Management**: Ensure the `GROQ_API_KEY` environment variable is set in production environments\n")
    
    print(f"Report generated: {report_file}")
    return str(report_file)

if __name__ == "__main__":
    report_file = generate_report()
    print(f"\nReport generated: {report_file}")
    
    # Open the report file if on a system that supports it
    try:
        if sys.platform == "darwin":  # macOS
            subprocess.run(["open", report_file])
        elif sys.platform == "win32":  # Windows
            os.startfile(report_file)
        elif sys.platform == "linux":  # Linux
            subprocess.run(["xdg-open", report_file])
    except Exception as e:
        print(f"Could not open report file: {e}") 