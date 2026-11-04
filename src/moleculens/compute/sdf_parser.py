"""Minimal SDF parser without RDKit dependency.

Parses V2000 MOL/SDF files; V3000 is explicitly unsupported.
"""

import math
from dataclasses import dataclass

ELEMENTS = [
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
]


@dataclass(frozen=True)
class Atom:
    """Represents an atom with its coordinates and element symbol."""

    symbol: str
    x: float
    y: float
    z: float


@dataclass
class SDFMolecule:
    """Represents a molecule parsed from SDF format."""

    name: str
    atoms: list[Atom]
    comment: str = ""

    @property
    def num_atoms(self) -> int:
        return len(self.atoms)

    def to_psi4_geometry(self, charge: int = 0, multiplicity: int = 1) -> str:
        """Convert to Psi4 geometry string with proper settings.

        Returns geometry block with:
        - units angstrom
        - no_com (no center of mass shift)
        - no_reorient (no rotation)
        """
        lines = [
            f"{charge} {multiplicity}",
            "units angstrom",
            "no_com",
            "no_reorient",
            "",
        ]
        for atom in self.atoms:
            lines.append(f"{atom.symbol}  {atom.x:12.6f}  {atom.y:12.6f}  {atom.z:12.6f}")
        return "\n".join(lines)

    def geometry_hash(self) -> str:
        """Compute a stable hash of the geometry for caching."""
        import hashlib

        def _normalize_coordinate(value: float) -> float:
            rounded = round(value, 6)
            return 0.0 if rounded == -0.0 else rounded

        # Canonicalize atom order so equivalent SDFs reuse the same cache entry.
        data = []
        for symbol, x, y, z in sorted(
            (
                atom.symbol,
                _normalize_coordinate(atom.x),
                _normalize_coordinate(atom.y),
                _normalize_coordinate(atom.z),
            )
            for atom in self.atoms
        ):
            data.append(f"{symbol}:{x:.6f}:{y:.6f}:{z:.6f}")
        content = "|".join(data)
        return hashlib.sha256(content.encode()).hexdigest()[:16]


class SDFParseError(Exception):
    """Raised when SDF parsing fails."""


def parse_sdf(content: str) -> SDFMolecule:
    """Parse an SDF/MOL file content and extract molecule data.

    Supports V2000 format (standard MOL format).

    Args:
        content: The SDF file content as a string.

    Returns:
        SDFMolecule with parsed atom data.

    Raises:
        SDFParseError: If the file format is invalid or parsing fails.
    """
    if len(content.encode()) > 2_000_000:
        raise SDFParseError("SDF exceeds 2 MB limit")
    lines = content.replace("\r\n", "\n").replace("\r", "\n").split("\n")

    if len(lines) < 4:
        raise SDFParseError("SDF file too short - needs at least 4 header lines")

    # Line 1: Molecule name
    name = lines[0].strip()

    # Line 2: Program/timestamp (ignored)
    # Line 3: Comment
    comment = lines[2].strip() if len(lines) > 2 else ""

    # Line 4: Counts line - parse atom and bond counts
    counts_line = lines[3]

    # Check for V3000 format
    if "V3000" in counts_line:
        raise SDFParseError("V3000 format not yet supported - please use V2000 format")

    if "V2000" not in counts_line:
        raise SDFParseError("V2000 counts line required")

    # Parse V2000 counts line
    # Format: aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
    # aaa = number of atoms, bbb = number of bonds
    try:
        num_atoms = int(counts_line[0:3])
        num_bonds = int(counts_line[3:6])
    except (ValueError, IndexError) as e:
        raise SDFParseError(f"Invalid counts line: {counts_line!r}") from e

    if not 0 < num_atoms <= 999 or not 0 <= num_bonds <= 999:
        raise SDFParseError(f"Invalid atom count: {num_atoms}")

    if len(lines) < 4 + num_atoms + num_bonds:
        raise SDFParseError(
            f"File has {len(lines)} lines but needs at least {4 + num_atoms} for {num_atoms} atoms"
        )

    # Parse atom block (lines 5 to 4+num_atoms)
    atoms: list[Atom] = []
    for i in range(num_atoms):
        atom_line = lines[4 + i]
        try:
            atom = _parse_atom_line(atom_line)
            atoms.append(atom)
        except ValueError as e:
            raise SDFParseError(f"Invalid atom line {i + 1}: {atom_line!r}") from e

    for line in lines[4 + num_atoms : 4 + num_atoms + num_bonds]:
        try:
            if len(line) < 9 or not all(1 <= int(line[i : i + 3]) <= num_atoms for i in (0, 3)):
                raise ValueError("Invalid bond index")
        except ValueError as exc:
            raise SDFParseError("Invalid bond record") from exc
    if "M  END" not in [line.strip() for line in lines[4 + num_atoms + num_bonds :]]:
        raise SDFParseError("Missing M END")
    return SDFMolecule(name=name, atoms=atoms, comment=comment)


def _parse_atom_line(line: str) -> Atom:
    """Parse a single atom line from V2000 format.

    Format (fixed width):
    xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee

    Where:
    - x, y, z: Coordinates (10 chars each)
    - aaa: Atom symbol (3 chars)
    """
    if len(line) < 34:
        raise ValueError(f"Atom line too short: {len(line)} chars")
    x = float(line[0:10].strip())
    y = float(line[10:20].strip())
    z = float(line[20:30].strip())
    symbol = line[31:34].strip()

    if symbol not in ELEMENTS:
        raise ValueError("Unknown element symbol")
    if not all(math.isfinite(v) and abs(v) <= 1_000_000 for v in (x, y, z)):
        raise ValueError("Coordinates must be finite and within supported bounds")

    return Atom(symbol=symbol, x=x, y=y, z=z)
