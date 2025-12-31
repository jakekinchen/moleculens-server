"""Unit tests for SDF parser."""

import pytest

from moleculens.compute.sdf_parser import Atom, SDFMolecule, SDFParseError, parse_sdf


class TestSDFParser:
    """Tests for SDF file parsing."""

    def test_parse_water_sdf(self, water_sdf: str) -> None:
        """Test parsing water molecule SDF."""
        mol = parse_sdf(water_sdf)

        assert mol.name == "water"
        assert mol.num_atoms == 3

        # Check atoms
        assert mol.atoms[0].symbol == "O"
        assert mol.atoms[0].x == pytest.approx(0.0)
        assert mol.atoms[0].y == pytest.approx(0.0)
        assert mol.atoms[0].z == pytest.approx(0.0)

        assert mol.atoms[1].symbol == "H"
        assert mol.atoms[1].x == pytest.approx(0.7586)
        assert mol.atoms[1].z == pytest.approx(0.5043)

        assert mol.atoms[2].symbol == "H"
        assert mol.atoms[2].x == pytest.approx(-0.7586)

    def test_to_psi4_geometry(self, water_sdf: str) -> None:
        """Test conversion to Psi4 geometry string."""
        mol = parse_sdf(water_sdf)
        geom = mol.to_psi4_geometry()

        assert "units angstrom" in geom
        assert "no_com" in geom
        assert "no_reorient" in geom
        assert "O" in geom
        assert "H" in geom

    def test_geometry_hash_stable(self, water_sdf: str) -> None:
        """Test that geometry hash is stable."""
        mol1 = parse_sdf(water_sdf)
        mol2 = parse_sdf(water_sdf)

        assert mol1.geometry_hash() == mol2.geometry_hash()
        assert len(mol1.geometry_hash()) == 16  # Truncated SHA256

    def test_parse_error_empty(self) -> None:
        """Test error on empty input."""
        with pytest.raises(SDFParseError):
            parse_sdf("")

    def test_parse_error_invalid_counts(self) -> None:
        """Test error on invalid counts line."""
        invalid = "test\n\n\nxxx yyy zzz\n"
        with pytest.raises(SDFParseError):
            parse_sdf(invalid)

    def test_parse_methane(self) -> None:
        """Test parsing methane SDF."""
        methane_sdf = """methane
  test

  5  4  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6287    0.6287    0.6287 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6287   -0.6287    0.6287 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6287   -0.6287   -0.6287 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6287    0.6287   -0.6287 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END
$$"""
        mol = parse_sdf(methane_sdf)

        assert mol.name == "methane"
        assert mol.num_atoms == 5
        assert mol.atoms[0].symbol == "C"
        assert sum(1 for a in mol.atoms if a.symbol == "H") == 4


class TestAtom:
    """Tests for Atom dataclass."""

    def test_atom_creation(self) -> None:
        """Test creating an atom."""
        atom = Atom(symbol="O", x=1.0, y=2.0, z=3.0)
        assert atom.symbol == "O"
        assert atom.x == 1.0
        assert atom.y == 2.0
        assert atom.z == 3.0

    def test_atom_immutable(self) -> None:
        """Test that Atom is immutable."""
        from dataclasses import FrozenInstanceError

        atom = Atom(symbol="O", x=1.0, y=2.0, z=3.0)
        with pytest.raises(FrozenInstanceError):
            atom.x = 5.0  # type: ignore


class TestSDFMolecule:
    """Tests for SDFMolecule dataclass."""

    def test_num_atoms(self) -> None:
        """Test num_atoms property."""
        mol = SDFMolecule(
            name="test",
            atoms=[
                Atom("O", 0, 0, 0),
                Atom("H", 1, 0, 0),
            ],
        )
        assert mol.num_atoms == 2
