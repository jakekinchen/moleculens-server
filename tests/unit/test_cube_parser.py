"""Unit tests for cube file parser."""

from pathlib import Path
from textwrap import dedent

import numpy as np
import pytest

from moleculens.compute.cube_parser import (
    BOHR_TO_ANGSTROM,
    CubeParseError,
    parse_cube_file,
)


class TestCubeParser:
    """Tests for Gaussian cube file parsing."""

    @pytest.fixture
    def sample_cube_file(self, tmp_path: Path) -> Path:
        """Create a sample cube file for testing."""
        content = dedent("""\
            Test cube file
            Generated for testing
                1    0.000000    0.000000    0.000000
                3    0.500000    0.000000    0.000000
                3    0.000000    0.500000    0.000000
                3    0.000000    0.000000    0.500000
                8    8.000000    0.000000    0.000000    0.000000
             1.0  2.0  3.0  4.0  5.0  6.0  7.0  8.0  9.0
            10.0 11.0 12.0 13.0 14.0 15.0 16.0 17.0 18.0
            19.0 20.0 21.0 22.0 23.0 24.0 25.0 26.0 27.0
        """)

        cube_file = tmp_path / "test.cube"
        cube_file.write_text(content)
        return cube_file

    def test_parse_cube_file(self, sample_cube_file: Path) -> None:
        """Test parsing a simple cube file."""
        cube = parse_cube_file(sample_cube_file)

        assert cube.n_points == (3, 3, 3)
        assert cube.data.shape == (3, 3, 3)
        assert len(cube.atom_numbers) == 1
        assert cube.atom_numbers[0] == 8  # Oxygen

        # Check origin
        assert np.allclose(cube.origin, [0.0, 0.0, 0.0])

        # Check axis vectors
        assert np.allclose(cube.axes[0], [0.5, 0.0, 0.0])
        assert np.allclose(cube.axes[1], [0.0, 0.5, 0.0])
        assert np.allclose(cube.axes[2], [0.0, 0.0, 0.5])

        # Check data values
        assert cube.data[0, 0, 0] == pytest.approx(1.0)
        assert cube.data[2, 2, 2] == pytest.approx(27.0)

    def test_origin_angstrom(self, sample_cube_file: Path) -> None:
        """Test origin conversion to Angstrom."""
        cube = parse_cube_file(sample_cube_file)
        origin_ang = cube.origin_angstrom

        assert np.allclose(origin_ang, [0.0, 0.0, 0.0])

    def test_axes_angstrom(self, sample_cube_file: Path) -> None:
        """Test axis conversion to Angstrom."""
        cube = parse_cube_file(sample_cube_file)
        axes_ang = cube.axes_angstrom

        expected = cube.axes * BOHR_TO_ANGSTROM
        assert np.allclose(axes_ang, expected)

    def test_grid_to_angstrom(self, sample_cube_file: Path) -> None:
        """Test grid coordinate conversion."""
        cube = parse_cube_file(sample_cube_file)

        # Origin point
        pos = cube.grid_to_angstrom(0, 0, 0)
        assert np.allclose(pos, [0.0, 0.0, 0.0])

        # Point at (1, 1, 1)
        pos = cube.grid_to_angstrom(1, 1, 1)
        expected = (cube.origin + cube.axes[0] + cube.axes[1] + cube.axes[2]) * BOHR_TO_ANGSTROM
        assert np.allclose(pos, expected)

    def test_parse_error_not_found(self) -> None:
        """Test error on missing file."""
        with pytest.raises(CubeParseError, match="not found"):
            parse_cube_file("/nonexistent/file.cube")

    def test_parse_error_too_short(self, tmp_path: Path) -> None:
        """Test error on truncated file."""
        cube_file = tmp_path / "short.cube"
        cube_file.write_text("header\nonly\n")

        with pytest.raises(CubeParseError, match="too short"):
            parse_cube_file(cube_file)


class TestCubeData:
    """Tests for CubeData dataclass."""

    def test_bohr_to_angstrom_constant(self) -> None:
        """Test the Bohr to Angstrom conversion constant."""
        # 1 Bohr ≈ 0.529 Angstrom
        assert 0.52 < BOHR_TO_ANGSTROM < 0.54
