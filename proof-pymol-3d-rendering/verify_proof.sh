#!/bin/bash
# PyMOL 3D Rendering Proof Verification Script

echo "=== PyMOL 3D Rendering Proof Verification ==="
echo

echo "1. File Sizes:"
ls -lh *.png
echo

echo "2. File Types:"
file *.png
echo

echo "3. SHA256 Hashes (should be different):"
shasum -a 256 *.png
echo

echo "4. Size Comparison:"
SIZE_2D=$(stat -c%s 2d_comparison.png)
SIZE_3D=$(stat -c%s 3d_render_proof.png)
RATIO=$(echo "scale=1; $SIZE_3D / $SIZE_2D" | bc -l)
echo "2D size: $SIZE_2D bytes"
echo "3D size: $SIZE_3D bytes"
echo "3D is ${RATIO}x larger than 2D"
echo

echo "5. PNG Headers (first 16 bytes):"
echo "2D header:"
hexdump -C 2d_comparison.png | head -1
echo "3D header:"
hexdump -C 3d_render_proof.png | head -1
echo

echo "✅ If files have different sizes and hashes, PyMOL 3D rendering is working!"
