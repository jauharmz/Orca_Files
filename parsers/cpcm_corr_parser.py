"""
ORCA CPCM Correction File Parser (.cpcm_corr)

Parses CPCM correction data including:
- Corrected dielectric energy
- Total C-PCM charge
- Corrected surface charges
"""

from dataclasses import dataclass, field


@dataclass
class CPCMCorrData:
    """Parsed CPCM correction file data."""
    corrected_energy: float = 0.0
    total_charge: float = 0.0
    corrected_charges: list[float] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            'corrected_energy': self.corrected_energy,
            'corrected_energy_kcal': self.corrected_energy * 627.509,
            'total_charge': self.total_charge,
            'num_charges': len(self.corrected_charges),
            'charge_stats': {
                'min': min(self.corrected_charges) if self.corrected_charges else 0,
                'max': max(self.corrected_charges) if self.corrected_charges else 0,
                'sum': sum(self.corrected_charges)
            } if self.corrected_charges else None
        }


def parse_cpcm_corr_file(filepath: str) -> CPCMCorrData:
    """Parse ORCA CPCM correction file."""
    with open(filepath, 'r') as f:
        content = f.read()
    return parse_cpcm_corr_content(content)


def parse_cpcm_corr_content(content: str) -> CPCMCorrData:
    """Parse CPCM correction content from string."""
    data = CPCMCorrData()
    lines = content.strip().split('\n')

    in_charges = False

    for line in lines:
        line = line.strip()

        if 'Corrected dielectric energy' in line:
            parts = line.split('=')
            if len(parts) >= 2:
                data.corrected_energy = float(parts[1].strip())

        elif 'Total C-PCM charge' in line:
            parts = line.split('=')
            if len(parts) >= 2:
                data.total_charge = float(parts[1].strip())

        elif 'C-PCM corrected charges' in line:
            in_charges = True

        elif in_charges and line:
            try:
                charge = float(line)
                data.corrected_charges.append(charge)
            except ValueError:
                pass

    return data


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        data = parse_cpcm_corr_file(sys.argv[1])
        print(f"Corrected energy: {data.corrected_energy:.6f} Eh")
        print(f"                : {data.corrected_energy * 627.509:.2f} kcal/mol")
        print(f"Total charge: {data.total_charge:.9f}")
        print(f"Surface charges: {len(data.corrected_charges)}")

        if data.corrected_charges:
            print(f"  Min: {min(data.corrected_charges):.6f}")
            print(f"  Max: {max(data.corrected_charges):.6f}")
            print(f"  Sum: {sum(data.corrected_charges):.9f}")
