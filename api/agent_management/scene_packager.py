"""
Scene Packager - Utility to package geometry and animation into a complete Three.js scene.
"""

from typing import Dict

from agent_management.models import (
    AnimationCode,
    FinalScenePackage,
    OrchestrationPlan,
    SceneScript,
)


class ScenePackager:
    """Takes the output from the geometry and animation agents and packages it
    into a complete Three.js scene."""

    @staticmethod
    def create_scene_package(
        script: SceneScript,
        orchestration_plan: OrchestrationPlan,
        object_geometries: Dict[str, Dict],
        animation_code: AnimationCode,
    ) -> FinalScenePackage:
        """Create a complete Three.js scene from the generated components.

        Args:
            script: The scene script with timecodes, descriptions, and captions
            orchestration_plan: The plan of objects needed for the scene
            object_geometries: Dictionary of objects with their generated geometry
            animation_code: Generated animation code

        Returns:
            FinalScenePackage: Complete scene package with HTML and JS
        """
        # Extract successful geometries
        successful_geometries = {
            name: data
            for name, data in object_geometries.items()
            if name != "_summary" and data.get("status") == "success"
        }

        # Count total elements
        total_elements = len(successful_geometries)

        # Extract all geometry code snippets
        all_geometry_code = "\n\n".join(
            [
                data.get("code", "").strip()
                for name, data in successful_geometries.items()
            ]
        )

        # Extract timecode markers
        timecode_markers = [keyframe.timecode for keyframe in animation_code.keyframes]

        # Create scene title
        title = orchestration_plan.scene_title

        # Create the JavaScript code
        js = ScenePackager._create_js(title, all_geometry_code, animation_code.code)
        minimal_js = ScenePackager._create_minimal_js(
            all_geometry_code, animation_code.code
        )

        # Create the HTML with embedded JavaScript
        html = ScenePackager._create_html(title, script, js)

        return FinalScenePackage(
            html=html,
            js=js,
            minimal_js=minimal_js,
            title=title,
            timecode_markers=timecode_markers,
            total_elements=total_elements,
        )

    @staticmethod
    def _create_html(title: str, script: SceneScript, js_code: str = None) -> str:
        """Create the complete HTML file with the Three.js scene.

        Args:
            title (str): Scene title
            script (SceneScript): Script with captions
            js_code (str, optional): JavaScript code to embed directly in the HTML (if None, will use a script link)
        Returns:
            str: Complete HTML file
        """
        # Create caption divs from script points
        caption_divs = []
        for point in script.content:
            caption_div = f'        <div class="caption" data-time="{point.timecode}" style="display: none;">{point.caption}</div>'
            caption_divs.append(caption_div)

        caption_html = "\n".join(caption_divs)

        # Determine whether to use embedded JS or external script link
        if js_code:
            script_tag = f"""    <!-- Embedded scene script -->
    <script>
{js_code}
    </script>"""
        else:
            script_tag = """    <!-- Main scene script -->
    <script src="/static/scene.js"></script>"""

        # Create the HTML template with corrected script paths
        html_template = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <style>
        body {{
            margin: 0;
            overflow: hidden;
            font-family: Arial, sans-serif;
            color: white;
        }}
        canvas {{
            width: 100%;
            height: 100%;
            display: block;
        }}
        .caption-container {{
            position: absolute;
            bottom: 30px;
            left: 0;
            width: 100%;
            text-align: center;
            z-index: 100;
            pointer-events: none;
        }}
        .caption {{
            display: none;
            background-color: rgba(0, 0, 0, 0.5);
            padding: 10px 20px;
            border-radius: 5px;
            margin: 0 auto;
            max-width: 80%;
            font-weight: bold;
            text-shadow: 1px 1px 2px rgba(0,0,0,0.8);
        }}
        .title {{
            position: absolute;
            top: 20px;
            left: 0;
            width: 100%;
            text-align: center;
            z-index: 100;
            font-size: 24px;
            pointer-events: none;
            text-shadow: 1px 1px 2px rgba(0,0,0,0.8);
        }}
        .controls {{
            position: absolute;
            bottom: 10px;
            right: 10px;
            z-index: 100;
            display: flex;
            gap: 10px;
        }}
        .controls button {{
            background-color: rgba(0, 0, 0, 0.5);
            border: 1px solid white;
            color: white;
            padding: 5px 10px;
            cursor: pointer;
            border-radius: 4px;
        }}
        .controls button:hover {{
            background-color: rgba(0, 0, 0, 0.7);
        }}
        .timeline {{
            position: absolute;
            bottom: 0;
            left: 0;
            height: 5px;
            background-color: rgba(255, 255, 255, 0.3);
            width: 100%;
        }}
        .progress {{
            height: 100%;
            width: 0%;
            background-color: rgba(255, 255, 255, 0.8);
        }}
    </style>
</head>
<body>
    <div class="title">{title}</div>
    <div class="caption-container">
{caption_html}
    </div>
    <div class="controls">
        <button id="rewind">⏪</button>
        <button id="play-pause">Pause</button>
        <button id="fast-forward">⏩</button>
        <button id="reset">Reset</button>
    </div>
    <div class="timeline">
        <div class="progress" id="progress-bar"></div>
    </div>

    <!-- Canvas element for the visualization -->
    <canvas id="scene-canvas" style="width: 100%; height: 100%; display: block;"></canvas>

    <!-- Main Three.js library -->
    <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/build/three.min.js"></script>
    <!-- OrbitControls add-on (corrected path) -->
    <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
    {script_tag}
    </body>
    </html>
    """
        return html_template

    @staticmethod
    def _create_js(title: str, geometry_code: str, animation_code: str) -> str:
        """Create the complete JavaScript file with the Three.js scene.

        Args:
            title: Scene title
            geometry_code: Combined geometry code from the geometry agent
            animation_code: Combined animation code from the animation agent

        Returns:
            str: Complete JavaScript file
        """
        # Create the JavaScript template header
        js_header = "// Scientific visualization: " + title + "\n"
        js_header += "// Auto-generated code\n\n"
        js_header += "// Global variables\n"
        js_header += "// Use existing globals if defined, otherwise create new ones\n"
        js_header += "var scene = window.scene || null;\n"
        js_header += "var camera = window.camera || null;\n"
        js_header += "var renderer = window.renderer || null;\n"
        js_header += "var controls = window.controls || null;\n"
        js_header += "var clock = window.clock || null;\n"
        js_header += "var isPlaying = typeof window.isPlaying !== 'undefined' ? window.isPlaying : true;\n"
        js_header += "// Set animation duration to 2 minutes\n"
        js_header += "var ANIMATION_DURATION = 120;\n\n"

        # Create the init function
        init_function = """
// Initialize the scene
function init() {
    // Create scene if it doesn't exist
    if (!scene) {
        scene = new THREE.Scene();
        scene.background = new THREE.Color(0x111122);
    }

    // Create camera if it doesn't exist
    if (!camera) {
        camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
        camera.position.z = 10;
    }

    // Look for an existing canvas or create renderer
    if (!renderer) {
        const canvas = document.getElementById('scene-canvas');
        if (canvas) {
            renderer = new THREE.WebGLRenderer({
                canvas: canvas,
                antialias: true
            });
            renderer.setSize(window.innerWidth, window.innerHeight);
            renderer.setPixelRatio(window.devicePixelRatio);
        } else {
            console.warn('No canvas element found with id "scene-canvas"');
            return; // Exit init if canvas not found
        }
    }

    // Create lighting
    const ambientLight = new THREE.AmbientLight(0x404040);
    scene.add(ambientLight);

    const directionalLight = new THREE.DirectionalLight(0xffffff, 1);
    directionalLight.position.set(1, 1, 1).normalize();
    scene.add(directionalLight);

    // Add orbit controls if they don't exist
    if (!controls && renderer) {
        controls = new THREE.OrbitControls(camera, renderer.domElement);
        controls.enableDamping = true;
        controls.dampingFactor = 0.25;
    }

    // Initialize clock if it doesn't exist
    if (!clock && typeof THREE !== 'undefined') {
        clock = new THREE.Clock();
    }

    // Create geometry
    createGeometry();

    // Handle window resize
    window.addEventListener('resize', onWindowResize);

    // Set up UI controls
    setupControls();
}
"""

        # Create the geometry function
        geometry_function = (
            """
// Create all geometry in the scene
function createGeometry() {
    // Geometry created by the GeometryAgent
"""
            + geometry_code
            + """
}
"""
        )

        # Create the animate function
        animate_function = (
            """
// Animation loop
function animate() {
    if (!isPlaying) {
        requestAnimationFrame(animate);
        renderer.render(scene, camera);
        controls.update();
        return;
    }

    requestAnimationFrame(animate);

    // Get the adjusted time with offset applied
    const adjustedTime = updateUI();

    // We'll make this time variable available to the animation code
    window.animationTime = adjustedTime;

    // Animation created by the AnimationAgent
"""
            + animation_code
            + """

    // Update controls
    controls.update();

    // Render the scene
    renderer.render(scene, camera);
}
"""
        )

        # Create the UI update function
        ui_function = """
// Update captions and UI based on current time
function updateUI() {
    // Add the time offset to the elapsed time
    const elapsedTime = clock.getElapsedTime() + (window.timeOffset || 0);

    // Constrain elapsed time to be between 0 and animation duration
    const constrainedTime = Math.max(0, Math.min(elapsedTime, ANIMATION_DURATION));

    // Update progress bar
    const progressBar = document.getElementById('progress-bar');
    const progress = constrainedTime / ANIMATION_DURATION;
    progressBar.style.width = progress * 100 + '%';

    // Update captions
    const captions = document.querySelectorAll('.caption');
    captions.forEach(caption => {
        const timeStr = caption.getAttribute('data-time');
        const [min, sec] = timeStr.split(':').map(Number);
        const timeInSeconds = min * 60 + sec;

        // Show caption if we're within 5 seconds of its timecode
        if (constrainedTime >= timeInSeconds && constrainedTime < timeInSeconds + 5) {
            caption.style.display = 'block';
        } else {
            caption.style.display = 'none';
        }
    });

    return constrainedTime; // Return the constrained time for use in animations
}
"""

        # Create the resize handler function
        resize_function = """
// Window resize handler
function onWindowResize() {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight);
}
"""

        # Create the controls setup function
        controls_function = """
// Set up UI controls
function setupControls() {
    // Store button references
    var playPauseButton = document.getElementById('play-pause');
    var resetButton = document.getElementById('reset');
    var rewindButton = document.getElementById('rewind');
    var fastForwardButton = document.getElementById('fast-forward');

    // Only set up event listeners if elements exist
    if (!playPauseButton || !resetButton || !rewindButton || !fastForwardButton) {
        console.warn('Some control buttons not found in the DOM. Control setup may be incomplete.');
    }

    let playbackSpeed = 1.0; // Normal speed

    // Initialize time offset if not already defined
    window.timeOffset = window.timeOffset || 0;

    // Set up play/pause button if it exists
    if (playPauseButton) {
        playPauseButton.addEventListener('click', () => {
            isPlaying = !isPlaying;
            playPauseButton.textContent = isPlaying ? 'Pause' : 'Play';
            if (isPlaying) {
                if (clock) clock.start();
            } else {
                if (clock) clock.stop();
            }
        });
    }

    // Set up reset button if it exists
    if (resetButton) {
        resetButton.addEventListener('click', () => {
            if (typeof THREE !== 'undefined') {
                clock = new THREE.Clock();
            }
            window.timeOffset = 0; // Reset the time offset
            isPlaying = true;
            if (playPauseButton) playPauseButton.textContent = 'Pause';
            playbackSpeed = 1.0; // Reset speed to normal
        });
    }

    // Set up rewind button if it exists
    if (rewindButton) {
        rewindButton.addEventListener('click', () => {
            // Decrease the time offset by 10 seconds (but don't go below negative animation duration)
            window.timeOffset = Math.max(window.timeOffset - 10, -ANIMATION_DURATION);

            // Ensure playing state
            isPlaying = true;
            if (playPauseButton) playPauseButton.textContent = 'Pause';
        });
    }

    // Set up fast-forward button if it exists
    if (fastForwardButton) {
        fastForwardButton.addEventListener('click', () => {
            // Increase the time offset by 10 seconds (but don't exceed animation duration)
            window.timeOffset = Math.min(window.timeOffset + 10, ANIMATION_DURATION);

            // Ensure playing state
            isPlaying = true;
            if (playPauseButton) playPauseButton.textContent = 'Pause';
        });
    }
}
"""

        # Create the initialization code
        init_code = """
// Initialize and start animation
init();
animate();
"""

        # Combine all parts into the complete JavaScript
        js_template = (
            js_header
            + init_function
            + geometry_function
            + animate_function
            + ui_function
            + resize_function
            + controls_function
            + init_code
        )

        return js_template

    @staticmethod
    def _create_minimal_js(geometry_code: str, animation_code: str) -> str:
        """Create minimal JavaScript code without boilerplate, for embedding in
        existing Three.js scenes.

        Args:
            geometry_code: Combined geometry code from the geometry agent
            animation_code: Combined animation code from the animation agent

        Returns:
            str: Minimal JavaScript code for embedding
        """
        # Create a minimal scene setup that can be injected into an existing canvas
        minimal_js = """// Function to create and run a Three.js scene in a given canvas element
function createScene(canvasElement) {
    // Scene setup
    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0x111122);

    // Camera setup
    const camera = new THREE.PerspectiveCamera(75, canvasElement.clientWidth / canvasElement.clientHeight, 0.1, 1000);
    camera.position.z = 10;

    // Renderer setup using the provided canvas
    const renderer = new THREE.WebGLRenderer({
        canvas: canvasElement,
        antialias: true
    });
    renderer.setSize(canvasElement.clientWidth, canvasElement.clientHeight);
    renderer.setPixelRatio(window.devicePixelRatio);

    // Lighting setup
    const ambientLight = new THREE.AmbientLight(0x404040);
    scene.add(ambientLight);

    const directionalLight = new THREE.DirectionalLight(0xffffff, 1);
    directionalLight.position.set(1, 1, 1).normalize();
    scene.add(directionalLight);

    // Controls setup
    const controls = new THREE.OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;
    controls.dampingFactor = 0.25;

    // Initialize clock and animation variables
    const clock = new THREE.Clock();
    let isPlaying = true;
    window.timeOffset = 0;

    // Create geometry objects
"""

        # Add the geometry code
        minimal_js += geometry_code + "\n\n"

        # Add the animation function with proper closure
        minimal_js += """    // Animation function
    function animate() {
        if (!isPlaying) {
            requestAnimationFrame(animate);
            renderer.render(scene, camera);
            controls.update();
            return;
        }

        requestAnimationFrame(animate);

        // Get the adjusted time with offset applied
        const elapsedTime = clock.getElapsedTime() + window.timeOffset;
        const ANIMATION_DURATION = 120; // 2 minutes in seconds
        const constrainedTime = Math.max(0, Math.min(elapsedTime, ANIMATION_DURATION));

        // Make this time variable available to the animation code
        window.animationTime = constrainedTime;

        // Animation code
"""

        # Add the animation code
        minimal_js += animation_code + "\n"

        # Add the render call and the function's closing braces
        minimal_js += """
        // Update controls and render
        controls.update();
        renderer.render(scene, camera);
    }

    // Add window resize handler
    function handleResize() {
        if (!canvasElement) return;

        camera.aspect = canvasElement.clientWidth / canvasElement.clientHeight;
        camera.updateProjectionMatrix();
        renderer.setSize(canvasElement.clientWidth, canvasElement.clientHeight);
    }

    // Add event listener for window resize
    window.addEventListener('resize', handleResize);

    // Start animation
    animate();

    // Return control functions for external use
    return {
        play: function() {
            isPlaying = true;
            clock.start();
        },
        pause: function() {
            isPlaying = false;
            clock.stop();
        },
        reset: function() {
            clock = new THREE.Clock();
            window.timeOffset = 0;
            isPlaying = true;
        },
        fastForward: function(seconds = 10) {
            const ANIMATION_DURATION = 120; // 2 minutes in seconds
            window.timeOffset = Math.min(window.timeOffset + seconds, ANIMATION_DURATION);
        },
        rewind: function(seconds = 10) {
            const ANIMATION_DURATION = 120; // 2 minutes in seconds
            window.timeOffset = Math.max(window.timeOffset - seconds, -ANIMATION_DURATION);
        },
        setTime: function(seconds) {
            window.timeOffset = clock.getElapsedTime() - seconds;
        },
        dispose: function() {
            window.removeEventListener('resize', handleResize);
            controls.dispose();
            renderer.dispose();
        }
    };
}

// Usage example:
// const canvas = document.getElementById('scene-canvas');
// const sceneController = createScene(canvas);
//
// To control the scene:
// sceneController.pause();
// sceneController.play();
// sceneController.fastForward();
// sceneController.rewind();
// sceneController.reset();
"""

        return minimal_js
