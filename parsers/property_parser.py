"""
ORCA Property File Parser (.property.txt)

Parses property files containing various computed properties.
"""

import re
from dataclasses import dataclass, field


@dataclass
class PropertyData:
    """Parsed property file data."""
    geometry_index: int = 0
    num_atoms: int = 0
    status: str = ""
    mulliken_charges: list[float] = field(default_factory=list)
    loewdin_charges: list[float] = field(default_factory=list)
    dipole_moment: list[float] = field(default_factory=list)
    orbital_energies: list[float] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            'geometry_index': self.geometry_index,
            'num_atoms': self.num_atoms,
            'status': self.status,
            'mulliken_charges': self.mulliken_charges,
            'loewdin_charges': self.loewdin_charges,
            'dipole_moment': self.dipole_moment,
            'orbital_energies': self.orbital_energies
        }


def parse_property_file(filepath: str) -> PropertyData:
    """Parse ORCA property file."""
    with open(filepath, 'r') as f:
        content = f.read()
    return parse_property_content(content)


def parse_property_content(content: str) -> PropertyData:
    """Parse property content from string."""
    data = PropertyData()

    # Status
    status_match = re.search(r'&Status\s+\[&Type "String"\]\s+"([^"]+)"', content)
    if status_match:
        data.status = status_match.group(1)

    # Geometry index
    geom_match = re.search(r'&GeometryIndex\s+(\d+)', content)
    if geom_match:
        data.geometry_index = int(geom_match.group(1))

    # Number of atoms
    natoms_match = re.search(r'&NAtoms\s+\[&Type "Integer"\]\s+(\d+)', content)
    if natoms_match:
        data.num_atoms = int(natoms_match.group(1))

    # Mulliken charges
    mulliken_section = re.search(
        r'\$SCF_Mulliken_Population_Analysis.*?&AtomicCharges.*?\n\s*\d+\s*\n(.*?)(?=\s*&|\$End)',
        content, re.DOTALL
    )
    if mulliken_section:
        for line in mulliken_section.group(1).strip().split('\n'):
            parts = line.split()
            if len(parts) >= 2:
                try:
                    data.mulliken_charges.append(float(parts[1]))
                except ValueError:
                    pass

    # Loewdin charges
    loewdin_section = re.search(
        r'\$SCF_Loewdin_Population_Analysis.*?&AtomicCharges.*?\n\s*\d+\s*\n(.*?)(?=\s*&|\$End)',
        content, re.DOTALL
    )
    if loewdin_section:
        for line in loewdin_section.group(1).strip().split('\n'):
            parts = line.split()
            if len(parts) >= 2:
                try:
                    data.loewdin_charges.append(float(parts[1]))
                except ValueError:
                    pass

    return data


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        data = parse_property_file(sys.argv[1])
        print(f"Status: {data.status}")
        print(f"Atoms: {data.num_atoms}")
        print(f"Mulliken charges: {len(data.mulliken_charges)}")
        print(f"Loewdin charges: {len(data.loewdin_charges)}")
