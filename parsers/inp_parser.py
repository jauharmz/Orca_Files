"""
ORCA Input File Parser (.inp)

Parses ORCA input files to extract calculation settings.
"""

import re
from dataclasses import dataclass, field


@dataclass
class InpData:
    """Parsed input file data."""
    keywords: list[str] = field(default_factory=list)
    method: str = ""
    basis_set: str = ""
    charge: int = 0
    multiplicity: int = 1
    job_types: list[str] = field(default_factory=list)
    maxcore: int = 0
    nprocs: int = 1
    coordinates: list[tuple[str, float, float, float]] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            'keywords': self.keywords,
            'method': self.method,
            'basis_set': self.basis_set,
            'charge': self.charge,
            'multiplicity': self.multiplicity,
            'job_types': self.job_types,
            'maxcore': self.maxcore,
            'nprocs': self.nprocs,
            'num_atoms': len(self.coordinates),
            'coordinates': [
                {'element': c[0], 'x': c[1], 'y': c[2], 'z': c[3]}
                for c in self.coordinates
            ]
        }


def parse_inp_file(filepath: str) -> InpData:
    """Parse ORCA input file."""
    with open(filepath, 'r') as f:
        content = f.read()
    return parse_inp_content(content)


def parse_inp_content(content: str) -> InpData:
    """Parse input content from string."""
    data = InpData()

    # Parse keyword line (starts with !)
    keyword_match = re.search(r'^!\s*(.+)$', content, re.MULTILINE)
    if keyword_match:
        keywords = keyword_match.group(1).split()
        data.keywords = keywords

        # Known methods
        methods = ['B3LYP', 'HF', 'MP2', 'CCSD', 'PBE', 'PBE0', 'M06', 'wB97X', 'BP86', 'TPSS']
        job_types = ['OPT', 'FREQ', 'SP', 'TD', 'TDDFT', 'NMR', 'NUMFREQ', 'RIJCOSX']

        for kw in keywords:
            kw_upper = kw.upper()
            # Check method
            for method in methods:
                if method in kw_upper:
                    data.method = kw
                    break
            # Check job type
            for jt in job_types:
                if jt in kw_upper:
                    data.job_types.append(jt)
            # Check basis set (contains numbers and characters like 6-31G)
            if re.match(r'\d+-\d+', kw) or 'def2' in kw.lower() or 'cc-pv' in kw.lower():
                data.basis_set = kw

    # Parse maxcore
    maxcore_match = re.search(r'%\s*maxcore\s+(\d+)', content, re.IGNORECASE)
    if maxcore_match:
        data.maxcore = int(maxcore_match.group(1))

    # Parse nprocs/PAL
    pal_match = re.search(r'PAL(\d+)', content)
    if pal_match:
        data.nprocs = int(pal_match.group(1))

    # Parse geometry block
    geom_match = re.search(r'\*\s*xyz\s+(-?\d+)\s+(\d+)\s*\n(.*?)\*', content, re.DOTALL)
    if geom_match:
        data.charge = int(geom_match.group(1))
        data.multiplicity = int(geom_match.group(2))

        for line in geom_match.group(3).strip().split('\n'):
            parts = line.split()
            if len(parts) >= 4:
                try:
                    coord = (parts[0], float(parts[1]), float(parts[2]), float(parts[3]))
                    data.coordinates.append(coord)
                except ValueError:
                    pass

    return data


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        data = parse_inp_file(sys.argv[1])
        print(f"Method: {data.method}")
        print(f"Basis: {data.basis_set}")
        print(f"Charge: {data.charge}, Mult: {data.multiplicity}")
        print(f"Job types: {', '.join(data.job_types)}")
        print(f"Atoms: {len(data.coordinates)}")
