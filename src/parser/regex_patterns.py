"""Shared regex patterns for ORCA parsing."""

import re

# Geometry patterns
CARTESIAN_BLOCK = re.compile(
    r"CARTESIAN COORDINATES \(ANGSTROEM\)\n-+\n(.*?)\n\n",
    re.DOTALL
)

INTERNAL_BLOCK = re.compile(
    r"INTERNAL COORDINATES \(ANGSTROEM\)\n-+\n(.*?)\n\n",
    re.DOTALL
)

GEOMETRY_FILE = re.compile(
    r"The coordinates will be read from file:\s*(\S+)"
)

XYZFILE_INPUT = re.compile(
    r"\*\s*xyzfile\s+(\d+)\s+(\d+)\s+(\S+)",
    re.IGNORECASE
)

# Energy patterns
GIBBS_ENERGY = re.compile(
    r"Final Gibbs free energy\s*\.{3,}\s+(-?\d+\.\d+)\s+Eh"
)

SINGLE_POINT_ENERGY = re.compile(
    r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)"
)

# Orbital patterns
ORBITAL_BLOCK = re.compile(
    r"ORBITAL ENERGIES\n-+\n\s*NO.*?\n(.*?)(?:\n\n|\Z)",
    re.DOTALL
)

SPIN_UP_BLOCK = re.compile(
    r"SPIN UP ORBITALS\n.*?\n(.*?)(?=\n\s*SPIN DOWN ORBITALS|\Z)",
    re.DOTALL
)

SPIN_DOWN_BLOCK = re.compile(
    r"SPIN DOWN ORBITALS\n.*?\n(.*?)(?:\n\n|\Z)",
    re.DOTALL
)

# Spectroscopy patterns
IR_SPECTRUM = re.compile(
    r"-{5,}\nIR SPECTRUM\n-+\n.*?\n-+\n(.*?)(?:\n\n|\Z)",
    re.DOTALL
)

RAMAN_SPECTRUM = re.compile(
    r"RAMAN SPECTRUM\s*-+\s*Mode.*?-+\s*(.*?)(?:\n\n|\Z)",
    re.DOTALL
)

VIBRATIONS = re.compile(
    r"VIBRATIONAL FREQUENCIES\n-+\n.*?\n\n(.*?)(?:\n\n|\Z)",
    re.DOTALL
)

NMR_SHIELDING = re.compile(
    r"CHEMICAL SHIELDING SUMMARY \(ppm\)\s*-+\s*Nucleus\s+Element\s+Isotropic\s+Anisotropy\s*-+\s*([\s\S]*?)(?=\n{2,}|\Z)",
    re.DOTALL
)

NMR_COUPLING = re.compile(
    r"SUMMARY OF ISOTROPIC COUPLING CONSTANTS\s*J\s*\(Hz\)\s*-+\s*([\s\S]*?)(?=\n{2,}|\Z)",
    re.DOTALL
)

# TD-DFT patterns
TDDFT_STATE = re.compile(
    r"STATE\s+(\d+):\s+E=\s+([-]?\d*\.?\d+[Ee]?[+-]?\d*)\s+au\s+([-]?\d*\.?\d+[Ee]?[+-]?\d*)\s+eV\s+([-]?\d*\.?\d+[Ee]?[+-]?\d*)\s+cm\*{1,2}-1([\s\S]+?)(?=STATE|\Z)"
)

TDDFT_TRANSITION = re.compile(
    r"(\d+)a\s*->\s*(\d+)a\s*:\s*([-]?\d*\.?\d+[Ee]?[+-]?\d*)\s*\(c=\s*([-]?\d*\.?\d+[Ee]?[+-]?\d*)\)"
)

# Dipole spectrum patterns
ELECTRIC_DIPOLE = re.compile(
    r"-+\n\s*ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS\s*\n-+\n.*?Transition.*?\n-+\n(.*?)(?:\n{2,}|\Z)",
    re.DOTALL | re.IGNORECASE
)

VELOCITY_DIPOLE = re.compile(
    r"-+\n\s*ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS\s*\n-+\n.*?Transition.*?\n-+\n(.*?)(?:\n{2,}|\Z)",
    re.DOTALL | re.IGNORECASE
)

# Mulliken patterns
MULLIKEN_CHARGES = re.compile(
    r"-+\n\s*MULLIKEN POPULATION ANALYSIS\s*\n-+\n\s*#\s*Atom\s*Element\s*Pop\s*Charge.*?\n\s*-+\n(.*?)(?:\n{2,}|\Z)",
    re.DOTALL | re.IGNORECASE
)

# Input section
INPUT_SECTION = re.compile(
    r"={10,}\s*(INPUT FILE|INPUT)\s*\n={10,}\n(.*?)\*{4}END OF INPUT\*{4}",
    re.DOTALL | re.IGNORECASE
)

# ESD flag
ESD_FLAG = re.compile(r"ESDFlag\s+(\w+)", re.IGNORECASE)
HESS_FLAG = re.compile(r"HESSFLAG\s+(\w+)", re.IGNORECASE)
