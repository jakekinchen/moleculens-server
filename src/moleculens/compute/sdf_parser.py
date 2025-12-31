"""Minimal SDF parser without RDKit dependency.

Parses V2000 and V3000 MOL/SDF files to extract atom coordinates and symbols.
"""

from dataclasses import dataclass


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

    def to_psi4_geometry(self) -> str:
        """Convert to Psi4 geometry string with proper settings.

        Returns geometry block with:
        - units angstrom
        - no_com (no center of mass shift)
        - no_reorient (no rotation)
        """
        lines = [
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

        # Round coordinates to avoid floating point variations
        data = []
        for atom in self.atoms:
            data.append(f"{atom.symbol}:{atom.x:.6f}:{atom.y:.6f}:{atom.z:.6f}")
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
    lines = content.strip().split("\n")

    if len(lines) < 4:
        raise SDFParseError("SDF file too short - needs at least 4 header lines")

    # Line 1: Molecule name
    name = lines[0].strip()

    # Line 2: Program/timestamp (ignored)
    # Line 3: Comment
    comment = lines[2].strip() if len(lines) > 2 else ""

    # Line 4: Counts line - parse atom and bond counts
    counts_line = lines[3].strip()

    # Check for V3000 format
    if "V3000" in content:
        raise SDFParseError("V3000 format not yet supported - please use V2000 format")

    # Parse V2000 counts line
    # Format: aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
    # aaa = number of atoms, bbb = number of bonds
    try:
        # Handle both fixed-width and space-separated formats
        parts = counts_line.split()
        # Use first field if space-separated, otherwise fixed-width format
        num_atoms = int(parts[0]) if len(parts) >= 2 else int(counts_line[0:3].strip())
    except (ValueError, IndexError) as e:
        raise SDFParseError(f"Invalid counts line: {counts_line!r}") from e

    if num_atoms <= 0:
        raise SDFParseError(f"Invalid atom count: {num_atoms}")

    if len(lines) < 4 + num_atoms:
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

    return SDFMolecule(name=name, atoms=atoms, comment=comment)


def _parse_atom_line(line: str) -> Atom:
    """Parse a single atom line from V2000 format.

    Format (fixed width):
    xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee

    Where:
    - x, y, z: Coordinates (10 chars each)
    - aaa: Atom symbol (3 chars)
    """
    # Try space-separated first (more common in practice)
    parts = line.split()

    if len(parts) >= 4:
        x = float(parts[0])
        y = float(parts[1])
        z = float(parts[2])
        symbol = parts[3].strip()
    else:
        # Fall back to fixed-width parsing
        if len(line) < 34:
            raise ValueError(f"Atom line too short: {len(line)} chars")
        x = float(line[0:10].strip())
        y = float(line[10:20].strip())
        z = float(line[20:30].strip())
        symbol = line[31:34].strip()

    # Clean up symbol (remove any extra characters)
    symbol = "".join(c for c in symbol if c.isalpha())

    if not symbol:
        raise ValueError("Empty atom symbol")

    # Capitalize correctly (first letter upper, rest lower)
    symbol = symbol[0].upper() + symbol[1:].lower() if len(symbol) > 1 else symbol.upper()

    return Atom(symbol=symbol, x=x, y=y, z=z)
