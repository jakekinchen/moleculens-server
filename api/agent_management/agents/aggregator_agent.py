"""
Aggregator Agent - merges geometry code, animation code, and caption code into a single embedded Three.js HTML.
"""

class AggregatorAgent:
    def __init__(self):
        pass

    def combine_into_html(self, geometry_code: str, animation_code: str, caption_code: str) -> str:
        """
        Creates a single HTML string with placeholders for geometry, animation, and caption code.
        """
        base_html = r"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Multi-Agent Three.js Demo</title>
    <style>
        body, html {
            margin: 0; 
            padding: 0; 
            overflow: hidden;
            background: #000;
        }
        #sceneCanvas {
            display: block; 
        }
        #caption {
            position: absolute;
            bottom: 20px;
            width: 100%;
            text-align: center;
            color: white;
            font-family: Arial, sans-serif;
            font-size: 1.3em;
            text-shadow: 1px 1px 2px rgba(0,0,0,0.5);
            pointer-events: none;
            z-index: 1000;
            padding: 10px;
        }
    </style>
    <script type="importmap">
        {
            "imports": {
                "three": "https://cdn.jsdelivr.net/npm/three@0.162.0/build/three.module.js",
                "three/addons/": "https://cdn.jsdelivr.net/npm/three@0.162.0/examples/jsm/"
            }
        }
    </script>
</head>
<body>
    <div id="caption"></div>
    <canvas id="sceneCanvas"></canvas>
    <script type="module">
        import * as THREE from 'three';

        // Basic Scene Setup
        const canvas = document.getElementById('sceneCanvas');
        const renderer = new THREE.WebGLRenderer({
            canvas,
            antialias: true
        });
        renderer.setSize(window.innerWidth, window.innerHeight);
        renderer.setPixelRatio(window.devicePixelRatio);
        
        const scene = new THREE.Scene();
        scene.background = new THREE.Color(0x000000);

        const camera = new THREE.PerspectiveCamera(
            75, window.innerWidth / window.innerHeight, 0.1, 1000
        );
        camera.position.z = 5;

        // Basic lighting
        const ambientLight = new THREE.AmbientLight(0xffffff, 0.3);
        scene.add(ambientLight);
        const pointLight = new THREE.PointLight(0xffffff, 1);
        pointLight.position.set(5, 5, 5);
        scene.add(pointLight);

        // Responsive resize
        window.addEventListener('resize', () => {
            camera.aspect = window.innerWidth / window.innerHeight;
            camera.updateProjectionMatrix();
            renderer.setSize(window.innerWidth, window.innerHeight);
        });

        // Main animate loop
        function animate() {
            requestAnimationFrame(animate);
            renderer.render(scene, camera);
        }
        animate();

        // [GEOMETRY_CODE]
        // [ANIMATION_CODE]
        // [CAPTION_CODE]
    </script>
</body>
</html>
"""
        # Insert each snippet in place of placeholders
        final_html = (base_html
                      .replace("// [GEOMETRY_CODE]", geometry_code)
                      .replace("// [ANIMATION_CODE]", animation_code)
                      .replace("// [CAPTION_CODE]", caption_code))
        return final_html

    def aggregate(self, geometry_code: str, animation_code: str, caption_code: str) -> str:
        """
        Combine them into a single final HTML.
        """
        return self.combine_into_html(geometry_code, animation_code, caption_code)