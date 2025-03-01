import base64
import json
import os
import tempfile
from typing import List, Dict, Tuple, Optional, Union, Literal, Any
import subprocess
import asyncio
import dotenv

# External libraries
import openai
from playwright.async_api import async_playwright

# Configuration
MAX_ITERATIONS = 5

class CodeQualityTool:
    def __init__(self, openai_api_key: str):
        self.openai_client = openai.Client(api_key=openai_api_key)
    
    async def process_code(
        self, 
        code: str, 
        language: Literal["javascript", "python"], 
        code_type: Literal["standalone", "component"],
        visual_outcome: str,
        timecodes: Optional[List[int]] = None
    ) -> Dict:
        """
        Main entry point for the code quality tool.
        
        Args:
            code: The source code to evaluate
            language: "javascript" or "python"
            code_type: "standalone" or "component"
            visual_outcome: Description of expected visual outcome
            timecodes: Optional list of millisecond timestamps for screenshots
            
        Returns:
            Dict with results including success status, feedback, and any screenshots
        """
        # First stage: Run and check for errors
        build_result = await self._check_build(code, language, code_type, timecodes)
        
        if not build_result["success"]:
            # Generate XML diff and instructions for fixing
            diff_prompt = self._generate_diff_prompt(code, build_result["error"])
            return {
                "success": False,
                "stage": "build",
                "error": build_result["error"],
                "code": code,
                "diff_prompt": diff_prompt
            }
        
        # Second stage: Check visual outcome
        visual_result = await self._check_visual_outcome(
            build_result["screenshots"], 
            visual_outcome
        )
        
        return {
            "success": True,
            "stage": "complete",
            "code": code,
            "screenshots": build_result["screenshots"],
            "visual_feedback": visual_result["feedback"]
        }
    
    async def submit_diff_and_check(
        self, 
        original_code: str,
        diff: str, 
        language: Literal["javascript", "python"], 
        code_type: Literal["standalone", "component"],
        visual_outcome: str,
        timecodes: Optional[List[int]] = None
    ) -> Dict:
        """
        Apply XML diff to code and rerun the check process.
        
        Args:
            original_code: Original code with errors
            diff: XML diff to apply
            language: "javascript" or "python"
            code_type: "standalone" or "component"
            visual_outcome: Description of expected visual outcome
            timecodes: Optional list of timestamps for screenshots
            
        Returns:
            Results from rerunning the process
        """
        # Apply the diff to the code
        updated_code = self._apply_xml_diff(original_code, diff)
        
        # Re-run the entire process with the updated code
        return await self.process_code(
            updated_code, 
            language, 
            code_type, 
            visual_outcome, 
            timecodes
        )
    
    async def _check_build(
        self, 
        code: str, 
        language: Literal["javascript", "python"], 
        code_type: Literal["standalone", "component"],
        timecodes: Optional[List[int]] = None
    ) -> Dict:
        """
        Check if the code builds/runs correctly and capture screenshots if specified.
        
        Returns:
            Dict with success flag, error message (if any), and screenshot paths
        """
        if language == "javascript":
            if code_type == "component":
                # For ThreeJS components, we need to add boilerplate
                return await self._run_threejs_component(code, timecodes)
            else:
                # For standalone HTML, run directly
                return await self._run_html_standalone(code, timecodes)
        elif language == "python":
            # Run Python code
            return await self._run_python_code(code, code_type, timecodes)
    
    async def _run_threejs_component(self, code: str, timecodes: Optional[List[int]] = None) -> Dict:
        """Run a ThreeJS component with appropriate boilerplate."""
        # Create boilerplate HTML with ThreeJS imports
        boilerplate = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <meta charset="utf-8">
            <title>ThreeJS Test</title>
            <style>
                body {{ margin: 0; overflow: hidden; }}
                canvas {{ width: 100%; height: 100%; display: block; }}
            </style>
        </head>
        <body>
            <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
            <script>
            // Create the standard Three.js scene, camera, and renderer
            const scene = new THREE.Scene();
            const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
            const renderer = new THREE.WebGLRenderer();
            renderer.setSize(window.innerWidth, window.innerHeight);
            document.body.appendChild(renderer.domElement);
            
            // Position camera
            camera.position.z = 5;
            
            // USER COMPONENT CODE STARTS HERE
            {code}
            // USER COMPONENT CODE ENDS HERE
            
            // Add the group to the scene if it exists and is not already added
            if (typeof group !== 'undefined' && group instanceof THREE.Group) {{
                scene.add(group);
            }}
            
            // Animation loop
            function animate() {{
                requestAnimationFrame(animate);
                renderer.render(scene, camera);
            }}
            animate();
            </script>
        </body>
        </html>
        """
        
        # Save to temporary file and run with Playwright
        return await self._run_html_in_playwright(boilerplate, timecodes)
    
    async def _run_html_standalone(self, code: str, timecodes: Optional[List[int]] = None) -> Dict:
        """Run standalone HTML directly."""
        return await self._run_html_in_playwright(code, timecodes)
    
    async def _run_html_in_playwright(self, html_code: str, timecodes: Optional[List[int]] = None) -> Dict:
        """
        Run HTML code in Playwright and capture screenshots.
        
        Returns:
            Dict with success flag, error message (if any), and screenshot paths
        """
        screenshots = []
        async with async_playwright() as p:
            try:
                browser = await p.chromium.launch(headless=True)
                page = await browser.new_page()
                
                # Set content and wait for load
                await page.set_content(html_code)
                await page.wait_for_load_state("networkidle")
                
                # Create temp directory for screenshots
                with tempfile.TemporaryDirectory() as temp_dir:
                    if not timecodes:
                        # Default: wait 1s and capture single frame
                        await page.wait_for_timeout(1000)
                        screenshot_path = os.path.join(temp_dir, "render_0000ms.png")
                        await page.screenshot(path=screenshot_path)
                        screenshots.append(screenshot_path)
                    else:
                        # Take screenshot at each specified timecode
                        for tc in sorted(timecodes):
                            await page.wait_for_timeout(tc if tc == timecodes[0] else tc - timecodes[timecodes.index(tc)-1])
                            screenshot_path = os.path.join(temp_dir, f"render_{tc}ms.png")
                            await page.screenshot(path=screenshot_path)
                            screenshots.append(screenshot_path)
                
                await browser.close()
                return {
                    "success": True,
                    "screenshots": screenshots
                }
            except Exception as e:
                if browser:
                    await browser.close()
                return {
                    "success": False,
                    "error": str(e),
                    "screenshots": screenshots
                }
    
    async def _run_python_code(self, code: str, code_type: str, timecodes: Optional[List[int]] = None) -> Dict:
        """
        Run Python code (either standalone or component).
        
        Returns:
            Dict with success flag, error message (if any), and screenshot paths
        """
        with tempfile.NamedTemporaryFile(suffix=".py", delete=False) as temp_file:
            temp_file_path = temp_file.name
            temp_file.write(code.encode('utf-8'))
        
        try:
            # Run the Python code
            result = subprocess.run(
                ["python", temp_file_path], 
                capture_output=True,
                text=True,
                timeout=30
            )
            
            if result.returncode != 0:
                # Error running the code
                return {
                    "success": False,
                    "error": result.stderr,
                    "screenshots": []
                }
            
            # Check if the code generated any output files/screenshots
            # This depends on how the Python code is expected to work
            screenshot_dir = os.path.dirname(temp_file_path)
            screenshots = [
                os.path.join(screenshot_dir, f) 
                for f in os.listdir(screenshot_dir) 
                if f.startswith("render_") and f.endswith(".png")
            ]
            
            return {
                "success": True,
                "screenshots": screenshots
            }
        except subprocess.TimeoutExpired:
            return {
                "success": False,
                "error": "Code execution timed out after 30 seconds",
                "screenshots": []
            }
        except Exception as e:
            return {
                "success": False,
                "error": str(e),
                "screenshots": []
            }
        finally:
            # Clean up temp file
            if os.path.exists(temp_file_path):
                os.unlink(temp_file_path)
    
    async def _check_visual_outcome(self, screenshot_paths: List[str], expected_outcome: str) -> Dict:
        """
        Use GPT-4o to verify if screenshots match the expected visual outcome.
        
        Args:
            screenshot_paths: List of paths to screenshots
            expected_outcome: Description of expected visual outcome
            
        Returns:
            Dict with feedback from GPT-4o
        """
        # Create a list to hold all content parts
        content_parts = []
        
        # Start with the text prompt
        content_parts.append({
            "type": "text", 
            "text": f"Evaluate if these rendered images match the expected visual outcome: '{expected_outcome}'"
        })
        
        # Add all screenshots to the content parts
        for path in screenshot_paths:
            image_data = self._encode_image_for_gpt(path)
            content_parts.append(image_data)
        
        # Query GPT-4o using the correct message format
        response = self.openai_client.chat.completions.create(
            model="gpt-4o",
            messages=[
                {
                    "role": "system", 
                    "content": "You are a vision expert that evaluates if rendered images match expected visual outcomes. Provide clear feedback about what matches and what doesn't."
                },
                {
                    "role": "user", 
                    "content": content_parts
                }
            ]
        )
        
        feedback = response.choices[0].message.content
        
        return {
            "feedback": feedback
        }
    
    def _generate_diff_prompt(self, code: str, error_message: str) -> str:
        """
        Generate a prompt instructing the agent to return an XML diff to fix the code.
        
        Args:
            code: The original code with errors
            error_message: The error message from the build process
            
        Returns:
            XML diff prompt
        """
        return f"""
        The following code has errors:
        
        ```
        {code}
        ```
        
        Error message:
        ```
        {error_message}
        ```
        
        Please return an XML diff to fix the code. Format your response as follows:
        
        <xml_diff>
        <!-- XML diff instructions here -->
        </xml_diff>
        
        The diff should use minimal changes to fix the error.
        """
    
    def _apply_xml_diff(self, original_code: str, diff: str) -> str:
        """
        Apply an XML diff to the original code.
        
        Args:
            original_code: The original code
            diff: XML diff instructions
            
        Returns:
            Updated code with diff applied
        """
        # In a real implementation, this would use a proper XML diff tool
        # Here's a simplified version for demonstration
        
        # Extract the content between <xml_diff> tags
        import re
        diff_content = re.search(r'<xml_diff>(.*?)</xml_diff>', diff, re.DOTALL)
        
        if not diff_content:
            return original_code  # No valid diff found
        
        # Parse instructions from the diff
        # This is a placeholder - a real implementation would parse the XML diff format
        # and apply specific changes to the original code
        
        # For demonstration, let's assume the diff contains replace instructions
        # like: <replace line="23">new code</replace>
        lines = original_code.split('\n')
        replace_pattern = r'<replace line="(\d+)">(.*?)</replace>'
        
        for match in re.finditer(replace_pattern, diff_content.group(1), re.DOTALL):
            line_num = int(match.group(1)) - 1  # 0-indexed
            new_line = match.group(2).strip()
            
            if 0 <= line_num < len(lines):
                lines[line_num] = new_line
        
        return '\n'.join(lines)
    
    def _encode_image_for_gpt(self, image_path: str) -> Dict:
        """
        Encode an image for use in GPT-4o API calls.
        
        Args:
            image_path: Path to the image file
            
        Returns:
            Dict with image content for API
        """
        with open(image_path, "rb") as img_file:
            b64_data = base64.b64encode(img_file.read()).decode('utf-8')
        
        return {
            "type": "image_url", 
            "image_url": {"url": f"data:image/png;base64,{b64_data}"}
        }

# Example usage:
async def main():
    tool = CodeQualityTool(openai_api_key=os.getenv("OPENAI_API_KEY"))
    
    # Example: Process ThreeJS component
    js_code = """
    const group = new THREE.Group();
    const geometry = new THREE.BoxGeometry(1, 1, 1);
    const material = new THREE.MeshBasicMaterial({color: 0xff0000});
    const cube = new THREE.Mesh(geometry, material);
    group.add(cube);
    
    // Animation
    function animate() {
        requestAnimationFrame(animate);
        cube.rotation.x += 0.01;
        cube.rotation.y += 0.01;
    }
    animate();
    """
    
    result = await tool.process_code(
        code=js_code,
        language="javascript",
        code_type="component",
        visual_outcome="A red spinning cube",
        timecodes=[1000, 2000, 3000]  # Screenshots at 1s, 2s, 3s
    )
    
    if not result["success"]:
        print(f"Build error: {result['error']}")
        print(f"Diff prompt: {result['diff_prompt']}")
        
        # Simulate agent providing a diff
        diff = """
        <xml_diff>
        <replace line="3">const material = new THREE.MeshBasicMaterial({color: 0xff0000});</replace>
        </xml_diff>
        """
        
        # Apply diff and check again
        fixed_result = await tool.submit_diff_and_check(
            original_code=js_code,
            diff=diff,
            language="javascript",
            code_type="component",
            visual_outcome="A red spinning cube",
            timecodes=[1000, 2000, 3000]
        )
        
        print(f"After fix: {fixed_result['visual_feedback']}")
    else:
        print(f"Visual feedback: {result['visual_feedback']}")

if __name__ == "__main__":
    asyncio.run(main())