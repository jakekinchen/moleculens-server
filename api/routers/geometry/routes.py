from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

router = APIRouter(
    prefix="/geometry",
    tags=["Prompt"],
    responses={404: {"description": "Not found"}},
)

class ReturnGeometry(BaseModel):
    jsx: str

@router.get("/cube", response_model=ReturnGeometry)
def get_cube_geometry():
    cube_jsx = """
    <Canvas camera={{ position: [0, 0, 5], fov: 75 }} style={{ width: '100%', height: '100%' }}>
      <color attach="background" args={["#111"]} />
      <OrbitControls />
      <ambientLight intensity={0.5} />
      <pointLight position={[10, 10, 10]} />
      <mesh>
        <boxGeometry />
        <meshStandardMaterial color="#3b82f6" />
      </mesh>
    </Canvas>
    """
    return {"jsx": cube_jsx}
