"""Native charged xTB translation test independent of HTTP serialization."""

import json

import numpy as np

from moleculens.compute.charges import _calculate_xtb_subprocess
from moleculens.compute.sdf_parser import Atom, SDFMolecule

xyz = np.array(
    [[0, 0, 0], [0.6, 0.6, 0.6], [-0.6, -0.6, 0.6], [0.6, -0.6, -0.6], [-0.6, 0.6, -0.6]]
)
shift = np.array([3.0, -2.0, 1.0])


def molecule(points):
    return SDFMolecule(
        "ammonium", [Atom(s, *p) for s, p in zip(["N", "H", "H", "H", "H"], points, strict=True)]
    )


a = _calculate_xtb_subprocess(molecule(xyz), 1, 1)
b = _calculate_xtb_subprocess(molecule(xyz + shift), 1, 1)
np.testing.assert_allclose(a.charges.sum(), 1, atol=1e-5)
np.testing.assert_allclose(b.charges.sum(), 1, atol=1e-5)
np.testing.assert_allclose(b.dipole_origin - a.dipole_origin, shift, atol=1e-8)
np.testing.assert_allclose(a.dipole_vector, b.dipole_vector, atol=1e-5)
print(
    json.dumps(
        {
            "engine": "native xTB 6.7.1 CLI",
            "charge": 1,
            "multiplicity": 1,
            "chargeSum": float(a.charges.sum()),
            "translatedChargeSum": float(b.charges.sum()),
            "dipoleDebye": a.dipole_vector.tolist(),
            "translatedDipoleDebye": b.dipole_vector.tolist(),
            "originShiftAngstrom": (b.dipole_origin - a.dipole_origin).tolist(),
            "passed": True,
        },
        indent=2,
    )
)
