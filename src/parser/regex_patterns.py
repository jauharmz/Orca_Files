"""
Regex patterns for ORCA parsing - MATCHED TO ORIGINAL orca_praser.py
"""

import re

# Numeral regex
FLOAT_RE = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"

# Geometry patterns
RE_CARTESIAN_BLOCKS = re.compile(r"CARTESIAN COORDINATES \(ANGSTROEM\)\n-+\n(.*?)(?:\n\n|\Z)", flags=re.S)
RE_INTERNAL_BLOCKS = re.compile(r"INTERNAL COORDINATES \(ANGSTROEM\)\n-+\n(.*?)(?:\n\n|\Z)", flags=re.S)

# Input section
RE_INPUT_FILE = re.compile(r"={10,}\s+(INPUT FILE|INPUT)\s*\n={10,}\n(.*?)\*{4}END OF INPUT\*{4}", flags=re.S | re.I)

# Optimized parameters
RE_OPT_PARAMS = re.compile(
    r"THE OPTIMIZATION HAS CONVERGED.*?Optimized Parameters:.*?-+\n(.*?)(?:\n{2,}|\Z)",
    flags=re.S | re.I
)

# Energy patterns
RE_FINAL_GIBBS = re.compile(r"Final Gibbs free energy\s*[:.]*\s*(-?\d+\.\d+)\s+Eh", flags=re.I)
RE_FINAL_SPE = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", flags=re.I)

# Orbital patterns
RE_ORBITAL_BLOCK = re.compile(r"ORBITAL ENERGIES\n-+\n\s*NO.*?\n(.*?)(?:\n\n|\Z)", flags=re.S)
RE_SPIN_UP = re.compile(r"SPIN UP ORBITALS\n.*?\n(.*?)(?=\n\s*SPIN DOWN ORBITALS|\Z)", flags=re.S)
RE_SPIN_DOWN = re.compile(r"SPIN DOWN ORBITALS\n.*?\n(.*?)(?:\n\n|\Z)", flags=re.S)

# Vibrations
RE_VIBS = re.compile(r"VIBRATIONAL FREQUENCIES\n-+\n\n(.*?)(?:\n{2,}|\Z)", flags=re.S)

# IR spectrum
RE_IR = re.compile(r"-{5,}\nIR SPECTRUM\n-+\n.*?\n-+\n(.*?)(?:\n\n|\Z)", flags=re.S)

# Raman spectrum
RE_RAMAN = re.compile(r"RAMAN SPECTRUM\s*-+\s*Mode.*?-+\s*(.*?)(?:\n\n|\Z)", flags=re.S)

# Mulliken
RE_MULLIKEN = re.compile(
    r"-{5,}\s*MULLIKEN POPULATION ANALYSIS\s*-+\s*\d+.*?-+\n(.*?)(?:\n{2,}|\Z)",
    flags=re.S | re.I
)

# NMR
RE_NMR_SHIELD = re.compile(
    r"CHEMICAL SHIELDING SUMMARY.*?Nucleus.*?Element.*?Isotropic.*?Anisotropy.*?-+\n(.*?)(?:\n{2,}|\Z)",
    flags=re.S | re.I
)
RE_NMR_COUPLING = re.compile(
    r"SUMMARY OF ISOTROPIC COUPLING CONSTANTS.*?J.*?Hz.*?-+\n(.*?)(?:\n{2,}|\Z)",
    flags=re.S | re.I
)

# TD-DFT
RE_TDDFT_STATE = re.compile(
    r"STATE\s+(\d+):\s+E=\s*([-]?\d*\.?\d+[Ee]?[+-]?\d*)\s+au\s+"
    r"([-]?\d*\.?\d+[Ee]?[+-]?\d*)\s+eV\s+([-]?\d*\.?\d+[Ee]?[+-]?\d*)\s+cm\*{1,2}-1([\s\S]+?)(?=STATE|\Z)"
)
RE_TD_TRANS = re.compile(r"(\d+)a\s*->\s*(\d+)a\s*:\s*([-]?\d*\.?\d+[Ee]?[+-]?\d*)\s*\(c=\s*([-]?\d*\.?\d+[Ee]?[+-]?\d*)\)")

# ESD flags
RE_ESDFLAG = re.compile(r"ESDFlag\s+(\w+)", flags=re.I)
RE_HESSFLAG = re.compile(r"HESSFLAG\s+(\w+)", flags=re.I)
