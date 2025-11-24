"""
ORCA Property File Parser (.property.txt)

Parses property files containing:
- Calculation info
- SCF and DFT energies
- Mulliken/Loewdin charges
- Mayer bond orders
- Dipole moments
- CIS absorption spectrum (excitations)
- Solvation details
"""

import re
from dataclasses import dataclass, field


@dataclass
class Excitation:
    """Single electronic excitation."""
    energy_ev: float
    energy_cm: float
    wavelength_nm: float
    oscillator_strength: float


@dataclass
class SolvationInfo:
    """Solvation model information."""
    solvent: str = ""
    epsilon: float = 0.0
    refrac: float = 0.0
    energy: float = 0.0
    surface_area: float = 0.0
    num_points: int = 0


@dataclass
class DFTEnergy:
    """DFT energy components."""
    final_energy: float = 0.0
    exchange: float = 0.0
    correlation: float = 0.0
    xc_energy: float = 0.0
    vdw_correction: float = 0.0


@dataclass
class PropertyData:
    """Parsed property file data."""
    geometry_index: int = 0
    num_atoms: int = 0
    num_electrons: int = 0
    num_basis_functions: int = 0
    charge: int = 0
    multiplicity: int = 1
    status: str = ""
    scf_energy: float = 0.0
    dft_energy: DFTEnergy = None
    mulliken_charges: list[float] = field(default_factory=list)
    loewdin_charges: list[float] = field(default_factory=list)
    atomic_numbers: list[int] = field(default_factory=list)
    dipole_magnitude: float = 0.0
    dipole_total: list[float] = field(default_factory=list)  # [x, y, z]
    mayer_bond_orders: list[tuple[int, int, float]] = field(default_factory=list)
    excitations: list[Excitation] = field(default_factory=list)
    solvation: SolvationInfo = None

    def to_dict(self) -> dict:
        return {
            'geometry_index': self.geometry_index,
            'num_atoms': self.num_atoms,
            'num_electrons': self.num_electrons,
            'num_basis_functions': self.num_basis_functions,
            'charge': self.charge,
            'multiplicity': self.multiplicity,
            'status': self.status,
            'scf_energy': self.scf_energy,
            'dft_energy': {
                'final_energy': self.dft_energy.final_energy,
                'exchange': self.dft_energy.exchange,
                'correlation': self.dft_energy.correlation,
                'xc_energy': self.dft_energy.xc_energy,
                'vdw_correction': self.dft_energy.vdw_correction
            } if self.dft_energy else None,
            'mulliken_charges': self.mulliken_charges,
            'loewdin_charges': self.loewdin_charges,
            'dipole_magnitude': self.dipole_magnitude,
            'dipole_total': self.dipole_total,
            'mayer_bond_orders': [
                {'atom1': b[0], 'atom2': b[1], 'order': b[2]}
                for b in self.mayer_bond_orders
            ],
            'excitations': [
                {
                    'energy_ev': e.energy_ev,
                    'wavelength_nm': e.wavelength_nm,
                    'oscillator_strength': e.oscillator_strength
                }
                for e in self.excitations
            ],
            'solvation': {
                'solvent': self.solvation.solvent,
                'epsilon': self.solvation.epsilon,
                'energy': self.solvation.energy,
                'surface_area': self.solvation.surface_area
            } if self.solvation else None
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

    # Calculation info
    natoms_match = re.search(r'&NumOfAtoms\s+\[&Type "Integer"\]\s+(\d+)', content)
    if natoms_match:
        data.num_atoms = int(natoms_match.group(1))

    nelectrons_match = re.search(r'&NumOfElectrons\s+\[&Type "Integer"\]\s+(\d+)', content)
    if nelectrons_match:
        data.num_electrons = int(nelectrons_match.group(1))

    nbasis_match = re.search(r'&NumOfBasisFuncts\s+\[&Type "Integer"\]\s+(\d+)', content)
    if nbasis_match:
        data.num_basis_functions = int(nbasis_match.group(1))

    charge_match = re.search(r'\$Calculation_Info.*?&Charge\s+\[&Type "Integer"\]\s+(-?\d+)', content, re.DOTALL)
    if charge_match:
        data.charge = int(charge_match.group(1))

    mult_match = re.search(r'\$Calculation_Info.*?&Mult\s+\[&Type "Integer"\]\s+(\d+)', content, re.DOTALL)
    if mult_match:
        data.multiplicity = int(mult_match.group(1))

    # SCF Energy
    scf_section = re.search(
        r'\$SCF_Energy.*?&totalEnergy.*?\n\s*\d+\s+\d+\s*\n\s*(-?\d+\.?\d*)',
        content, re.DOTALL
    )
    if scf_section:
        data.scf_energy = float(scf_section.group(1))

    # DFT Energy components
    dft_section = re.search(r'\$DFT_Energy(.*?)(?=\$End|\$\w)', content, re.DOTALL)
    if dft_section:
        dft = DFTEnergy()
        dft_content = dft_section.group(1)

        final_match = re.search(r'&finalEn.*?(-?\d+\.?\d*)', dft_content)
        if final_match:
            dft.final_energy = float(final_match.group(1))

        exch_match = re.search(r'&eExchange.*?(-?\d+\.?\d*)', dft_content)
        if exch_match:
            dft.exchange = float(exch_match.group(1))

        corr_match = re.search(r'&eCorr.*?(-?\d+\.?\d*)', dft_content)
        if corr_match:
            dft.correlation = float(corr_match.group(1))

        xc_match = re.search(r'&eXC.*?(-?\d+\.?\d*)', dft_content)
        if xc_match:
            dft.xc_energy = float(xc_match.group(1))

        data.dft_energy = dft

    # VdW Correction
    vdw_match = re.search(r'\$VdW_Correction.*?&vdW.*?(-?\d+\.?\d*)', content, re.DOTALL)
    if vdw_match and data.dft_energy:
        data.dft_energy.vdw_correction = float(vdw_match.group(1))

    # Mulliken charges
    mulliken_section = re.search(
        r'\$SCF_Mulliken_Population_Analysis.*?&AtomicCharges.*?&Dim\s+\(\s*\d+\s*,\s*\d+\s*\)\s*(.*?)(?=\s*&|\$End)',
        content, re.DOTALL
    )
    if mulliken_section:
        for line in mulliken_section.group(1).strip().split('\n'):
            parts = line.split()
            for p in parts:
                try:
                    data.mulliken_charges.append(float(p))
                except ValueError:
                    pass

    # Loewdin charges
    loewdin_section = re.search(
        r'\$SCF_Loewdin_Population_Analysis.*?&AtomicCharges.*?&Dim\s+\(\s*\d+\s*,\s*\d+\s*\)\s*(.*?)(?=\s*&|\$End)',
        content, re.DOTALL
    )
    if loewdin_section:
        for line in loewdin_section.group(1).strip().split('\n'):
            parts = line.split()
            for p in parts:
                try:
                    data.loewdin_charges.append(float(p))
                except ValueError:
                    pass

    # Dipole moment
    dipole_section = re.search(
        r'\$SCF_Dipole_Moment.*?&dipoleMagnitude.*?(-?\d+\.?\d*).*?'
        r'&dipoleTotal.*?&Dim\s+\(\s*\d+\s*,\s*\d+\s*\)\s*(.*?)(?=\s*&|\$End)',
        content, re.DOTALL
    )
    if dipole_section:
        data.dipole_magnitude = float(dipole_section.group(1))
        dipole_lines = dipole_section.group(2).strip().split('\n')
        for line in dipole_lines:
            parts = line.split()
            for p in parts:
                try:
                    data.dipole_total.append(float(p))
                except ValueError:
                    pass

    # Mayer bond orders
    mayer_section = re.search(
        r'\$SCF_Mayer_Population_Analysis.*?&BondOrders.*?&Dim\s+\(\s*(\d+)\s*,.*?\)\s*(.*?)(?=\s*&ATNO)',
        content, re.DOTALL
    )
    if mayer_section:
        num_bonds = int(mayer_section.group(1))
        bond_values = []
        for line in mayer_section.group(2).strip().split('\n'):
            parts = line.split()
            for p in parts:
                try:
                    bond_values.append(float(p))
                except ValueError:
                    pass

        # Also need to get the atom pairs from &components
        components_match = re.search(
            r'\$SCF_Mayer_Population_Analysis.*?&components.*?&Dim\s+\(\s*\d+\s*,\s*\d+\s*\)\s*(.*?)(?=\s*&BondOrders)',
            content, re.DOTALL
        )
        if components_match:
            comp_values = []
            for line in components_match.group(1).strip().split('\n'):
                parts = line.split()
                for p in parts:
                    try:
                        comp_values.append(int(p))
                    except ValueError:
                        pass

            # Components are in groups of 4: atom1, atom2, Z1, Z2
            for i in range(min(len(bond_values), len(comp_values) // 4)):
                atom1 = comp_values[i * 4]
                atom2 = comp_values[i * 4 + 1]
                order = bond_values[i]
                data.mayer_bond_orders.append((atom1, atom2, order))

    # CIS Absorption Spectrum
    cis_section = re.search(
        r'\$CIS_Absorption_Spectrum.*?&NTrans\s+\[&Type "Integer"\]\s+(\d+).*?'
        r'&ExcitationEnergies.*?&Dim\s+\(\s*\d+\s*,\s*\d+\s*\)\s*(.*?)(?=\s*&|\$End)',
        content, re.DOTALL
    )
    if cis_section:
        num_trans = int(cis_section.group(1))
        values = []
        for line in cis_section.group(2).strip().split('\n'):
            parts = line.split()
            for p in parts:
                try:
                    values.append(float(p))
                except ValueError:
                    pass

        # Each excitation has 11 columns
        cols_per_exc = 11
        for i in range(num_trans):
            if (i + 1) * cols_per_exc <= len(values):
                idx = i * cols_per_exc
                exc = Excitation(
                    energy_ev=values[idx],
                    energy_cm=values[idx + 1],
                    wavelength_nm=values[idx + 2],
                    oscillator_strength=values[idx + 3]
                )
                data.excitations.append(exc)

    # Solvation details
    solv_section = re.search(r'\$Solvation_Details(.*?)(?=\$End|\$\w)', content, re.DOTALL)
    if solv_section:
        solv = SolvationInfo()
        solv_content = solv_section.group(1)

        solvent_match = re.search(r'&Solvent.*?"([^"]+)"', solv_content)
        if solvent_match:
            solv.solvent = solvent_match.group(1)

        eps_match = re.search(r'&Epsilon.*?(-?\d+\.?\d*)', solv_content)
        if eps_match:
            solv.epsilon = float(eps_match.group(1))

        refrac_match = re.search(r'&Refrac.*?(-?\d+\.?\d*)', solv_content)
        if refrac_match:
            solv.refrac = float(refrac_match.group(1))

        energy_match = re.search(r'&CPCMDielEnergy.*?(-?\d+\.?\d*)', solv_content)
        if energy_match:
            solv.energy = float(energy_match.group(1))

        area_match = re.search(r'&SurfaceArea.*?(-?\d+\.?\d*)', solv_content)
        if area_match:
            solv.surface_area = float(area_match.group(1))

        points_match = re.search(r'&NPoints.*?(\d+)', solv_content)
        if points_match:
            solv.num_points = int(points_match.group(1))

        data.solvation = solv

    return data


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        data = parse_property_file(sys.argv[1])
        print(f"Status: {data.status}")
        print(f"Atoms: {data.num_atoms}")
        print(f"Electrons: {data.num_electrons}")
        print(f"Basis functions: {data.num_basis_functions}")
        print(f"Charge: {data.charge}, Mult: {data.multiplicity}")

        if data.scf_energy:
            print(f"SCF Energy: {data.scf_energy:.6f} Eh")

        if data.dft_energy:
            print(f"DFT Energy: {data.dft_energy.final_energy:.6f} Eh")
            print(f"  XC: {data.dft_energy.xc_energy:.6f}")
            print(f"  VdW: {data.dft_energy.vdw_correction:.6f}")

        if data.dipole_magnitude:
            print(f"Dipole: {data.dipole_magnitude:.4f} a.u.")

        print(f"Mulliken charges: {len(data.mulliken_charges)}")
        print(f"Loewdin charges: {len(data.loewdin_charges)}")
        print(f"Bond orders: {len(data.mayer_bond_orders)}")

        if data.excitations:
            print(f"\nExcitations: {len(data.excitations)}")
            for i, e in enumerate(data.excitations[:3]):
                print(f"  {i+1}. {e.wavelength_nm:.1f} nm, f={e.oscillator_strength:.4f}")

        if data.solvation:
            print(f"\nSolvation: {data.solvation.solvent}")
            print(f"  Epsilon: {data.solvation.epsilon:.2f}")
            print(f"  Energy: {data.solvation.energy:.6f} Eh")
