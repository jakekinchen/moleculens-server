/* eslint-disable */
import React, { useState, useMemo, useEffect, useRef } from 'react';
import { Canvas, useThree, useFrame } from '@react-three/fiber';
import { OrbitControls } from '@react-three/drei';
import { ArrowsPointingInIcon, ArrowsPointingOutIcon } from '@heroicons/react/24/outline';
import * as THREE from 'three';
import { LoadingFacts } from './LoadingFacts';
import { OrbitControls as ThreeOrbitControls } from 'three/examples/jsm/controls/OrbitControls';

// CSS styles for atom labels
const atomLabelStyles = `
  .atom-label {
    color: #fff;
    font-family: sans-serif;
    padding: 2px;
    background: rgba(0,0,0,0.6);
    border-radius: 3px;
    font-size: 12px;
    pointer-events: none;
    user-select: none;
  }
`;

declare global {
  interface Window {
    PDBLoader?: any;
    labelRenderer?: any;
    labelRendererResizeListener?: boolean;
    CSS2DRenderer?: any;
    CSS2DObject?: any;
    THREE?: any;
    createAtomLabel?: (position: any, text: string) => any;
  }
}

// Helper to auto-fit camera to entire scene on first render
const CameraController = () => {
  const { camera, scene } = useThree();

  useEffect(() => {
    requestAnimationFrame(() => {
      const objects = scene.children.filter(child => !(child instanceof THREE.Light));
      if (!objects.length) return;

      const box = new THREE.Box3();
      objects.forEach(object => box.expandByObject(object));

      const size = box.getSize(new THREE.Vector3());
      const center = box.getCenter(new THREE.Vector3());

      const diagonal = Math.sqrt(size.x * size.x + size.y * size.y + size.z * size.z);
      const scaleFactor = Math.max(1.2, Math.log10(diagonal) * 0.8);
      const distance = diagonal * scaleFactor;

      // Position camera at a nice angle
      const theta = Math.PI / 4; // 45 deg
      const phi = Math.PI / 6;   // 30 deg
      camera.position.set(
        center.x + distance * Math.sin(theta) * Math.cos(phi),
        center.y + distance * Math.sin(phi),
        center.z + distance * Math.cos(theta) * Math.cos(phi)
      );

      camera.lookAt(center);

      // If orbit controls are used, update its target
      const controls = camera.userData.controls;
      if (controls) {
        controls.target.copy(center);
      }

      camera.near = distance * 0.01;
      camera.far = distance * 10;
      camera.updateProjectionMatrix();
    });
  }, [camera, scene]);

  return (
    <OrbitControls 
      makeDefault 
      autoRotate 
      autoRotateSpeed={1.5}
      enableDamping
      dampingFactor={0.05}
      minDistance={1}
      maxDistance={1000}
    />
  );
};

interface VisualizationPanelProps {
  script?: string;
  html?: string;
  title?: string;
  isLoading?: boolean;
  isInteractive?: boolean;
}

interface DynamicSceneComponentProps {
  code: string;
}

const DynamicSceneComponent: React.FC<DynamicSceneComponentProps> = ({ code }) => {
  const { scene, camera, gl: renderer } = useThree();
  const controls = useRef<ThreeOrbitControls | null>(null);
  const labelRendererRef = useRef<any>(null);

  // Inject atom-label styling
  useEffect(() => {
    if (!document.getElementById('atom-label-styles')) {
      const style = document.createElement('style');
      style.id = 'atom-label-styles';
      style.textContent = atomLabelStyles;
      document.head.appendChild(style);
    }
    return () => {
      const styleElement = document.getElementById('atom-label-styles');
      if (styleElement) styleElement.remove();
    };
  }, []);

  // Wrap the Python-generated code so we can call it
  const createVisualizationWrapper = (pythonGeneratedCode: string) => {
    return function setupVisualization(
      THREE: any,
      scene: THREE.Scene,
      camera: THREE.Camera,
      controls: any,
      labelRenderer: any
    ) {
      if (!window.CSS2DObject || !window.CSS2DRenderer) {
        console.error('CSS2D classes not available');
      }
      if (labelRenderer && !window.labelRenderer) {
        window.labelRenderer = labelRenderer;
      }

      const config = {
        camera,
        controls,
        labelRenderer: window.labelRenderer || labelRenderer,
        enableAnnotations: true
      };

      const createMoleculeVisualization = new Function(
        'THREE',
        'scene',
        'camera',
        'controls',
        'config',
        `
        if (typeof CSS2DObject === 'undefined' && window.CSS2DObject) {
          const CSS2DObject = window.CSS2DObject;
        }
        if (!window.createAtomLabel) {
          window.createAtomLabel = function(pos, text) {
            if (!window.CSS2DObject) return null;
            const div = document.createElement('div');
            div.className = 'atom-label';
            div.textContent = text;
            const lbl = new window.CSS2DObject(div);
            lbl.position.copy(pos);
            return lbl;
          };
        }
        ${pythonGeneratedCode}
        `
      );

      try {
        return createMoleculeVisualization(THREE, scene, camera, controls, config);
      } catch (error) {
        console.error('Visualization error:', error);
        throw error;
      }
    };
  };

  // Render loop for CSS2D
  useFrame(() => {
    if (labelRendererRef.current && camera && scene) {
      labelRendererRef.current.render(scene, camera);
      if (window.labelRenderer && window.labelRenderer !== labelRendererRef.current) {
        window.labelRenderer.render(scene, camera);
      }
    }
  });

  // Load the user code and set up scene
  useEffect(() => {
    if (!code || !scene || !camera || !renderer || !controls.current) return;

    (async function setupScene() {
      try {
        const { PDBLoader } = await import('three/addons/loaders/PDBLoader.js');
        const { CSS2DRenderer, CSS2DObject } = await import('three/addons/renderers/CSS2DRenderer.js');

        window.PDBLoader = PDBLoader;
        window.CSS2DRenderer = CSS2DRenderer;
        window.CSS2DObject = CSS2DObject;
        window.THREE = THREE;

        // Attach a CSS2D renderer
        const container = document.querySelector('#container');
        if (!container) return;

        // Clear any old labelRenderer
        if (labelRendererRef.current) {
          try {
            container.removeChild(labelRendererRef.current.domElement);
          } catch (e) {
            console.warn('Error removing existing labelRenderer:', e);
          }
          labelRendererRef.current = null;
        }

        const labelRenderer = new CSS2DRenderer();
        const rect = container.getBoundingClientRect();
        labelRenderer.setSize(rect.width, rect.height);
        labelRenderer.domElement.style.position = 'absolute';
        labelRenderer.domElement.style.top = '0px';
        labelRenderer.domElement.style.left = '0px';
        labelRenderer.domElement.style.width = '100%';
        labelRenderer.domElement.style.height = '100%';
        labelRenderer.domElement.style.pointerEvents = 'none';
        labelRenderer.domElement.style.zIndex = '10';
        container.appendChild(labelRenderer.domElement);

        labelRendererRef.current = labelRenderer;
        window.labelRenderer = labelRenderer;

        const handleResize = () => {
          const rect = container.getBoundingClientRect();
          labelRenderer.setSize(rect.width, rect.height);
        };
        window.addEventListener('resize', handleResize);

        // Invoke the user visualization code
        const setupVisualization = createVisualizationWrapper(code);
        setupVisualization(THREE, scene, camera, controls.current, labelRenderer);

        return () => {
          window.removeEventListener('resize', handleResize);
        };
      } catch (error) {
        console.error('Scene setup error:', error);
      }
    })();

    // Cleanup on unmount
    return () => {
      scene.children.slice().forEach(child => {
        if (!(child instanceof THREE.Light)) {
          scene.remove(child);
        }
      });
      delete window.PDBLoader;
      delete window.CSS2DRenderer;
      delete window.CSS2DObject;
      delete window.THREE;
      delete window.createAtomLabel;

      const container = document.querySelector('#container');
      if (container && labelRendererRef.current) {
        try {
          container.removeChild(labelRendererRef.current.domElement);
        } catch (e) {
          console.warn('Error removing labelRenderer:', e);
        }
        labelRendererRef.current = null;
        window.labelRenderer = null;
      }
    };
  }, [code, scene, camera, renderer, controls]);

  return (
    <>
      <ambientLight intensity={0.4} />
      <directionalLight position={[1, 1, 1]} intensity={1} />
      <directionalLight position={[-1, -1, -1]} intensity={0.4} />
      <OrbitControls
        ref={controls}
        makeDefault
        autoRotate
        autoRotateSpeed={1.5}
        enableDamping
        dampingFactor={0.05}
        minDistance={1}
        maxDistance={1000}
      />
    </>
  );
};

export const VisualizationPanel: React.FC<VisualizationPanelProps> = ({
  script,
  html,
  title,
  isLoading = false,
  isInteractive = false
}) => {
  const [isExpanded, setIsExpanded] = useState(false);
  const [isTransitioning, setIsTransitioning] = useState(false);

  const handleExpand = () => {
    setIsTransitioning(true);
    setIsExpanded(!isExpanded);
    setTimeout(() => setIsTransitioning(false), 300);
  };

  const DynamicScene = useMemo(() => {
    if (!script) return null;
    return () => <DynamicSceneComponent code={script} />;
  }, [script]);

  return (
    <div className="absolute inset-0">
      <div
        className={`transition-all duration-300 bg-gray-800 rounded-lg shadow-lg border border-gray-700 overflow-hidden
          ${isExpanded ? 'fixed inset-0 z-50 m-0' : 'absolute inset-0 m-2'}`}
      >
        {/* Expand/Collapse button */}
        <button
          onClick={handleExpand}
          className="absolute top-2 right-2 z-10 p-2 text-gray-400 hover:text-white transition-colors"
          aria-label={isExpanded ? 'Collapse visualization' : 'Expand visualization'}
        >
          {isExpanded ? (
            <ArrowsPointingInIcon className="w-6 h-6" />
          ) : (
            <ArrowsPointingOutIcon className="w-6 h-6" />
          )}
        </button>

        {/* Optional loading overlay */}
        {isLoading && (
          <div className="absolute inset-0 flex items-center justify-center bg-gray-800 bg-opacity-90 z-20">
            <LoadingFacts isVisible={isLoading} showFacts />
          </div>
        )}

        {/* Optional title (e.g., for non-interactive usage) */}
        {title && (
          <div className="absolute top-2 left-2 z-10 bg-gray-800/80 px-3 py-1.5 rounded-lg">
            <h3 className="text-sm font-medium text-white truncate max-w-[80%]">
              {title}
            </h3>
          </div>
        )}

        <div
          id="container"
          className="w-full h-full relative overflow-hidden"
          style={{
            boxSizing: 'border-box',
            padding: 0,
            margin: 0,
            transform: 'none',
            willChange: 'transform'
          }}
        >
          <Canvas
            style={{
              position: 'absolute',
              top: 0,
              left: 0,
              width: '100%',
              height: '100%',
              transform: 'none',
              willChange: 'transform'
            }}
          >
            {/* Camera auto-fit and controls */}
            <CameraController />
            {DynamicScene && <DynamicScene />}
          </Canvas>
        </div>
      </div>
    </div>
  );
};