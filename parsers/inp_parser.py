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
    solvent: str = ""
    dispersion: str = ""
    scf_convergence: str = ""
    grid: str = ""
    ri_approximation: str = ""
    coordinates: list[tuple[str, float, float, float]] = field(default_factory=list)
    blocks: dict[str, str] = field(default_factory=dict)  # % blocks

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
            'solvent': self.solvent,
            'dispersion': self.dispersion,
            'scf_convergence': self.scf_convergence,
            'grid': self.grid,
            'ri_approximation': self.ri_approximation,
            'num_atoms': len(self.coordinates),
            'coordinates': [
                {'element': c[0], 'x': c[1], 'y': c[2], 'z': c[3]}
                for c in self.coordinates
            ],
            'blocks': self.blocks
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
        job_types = ['OPT', 'FREQ', 'SP', 'TD', 'TDDFT', 'NMR', 'NUMFREQ']
        scf_settings = ['VERYTIGHTSCF', 'TIGHTSCF', 'NORMALSCF', 'LOOSESCF']
        grids = ['DEFGRID1', 'DEFGRID2', 'DEFGRID3', 'GRID4', 'GRID5', 'GRID6', 'GRID7']
        ri_methods = ['RIJCOSX', 'RIJK', 'RIJONX', 'NORI']
        dispersion = ['D3BJ', 'D3', 'D4', 'D3ZERO']

        for kw in keywords:
            kw_upper = kw.upper()

            # Check method
            for method in methods:
                if method in kw_upper:
                    data.method = kw
                    break

            # Check job type
            for jt in job_types:
                if jt == kw_upper:
                    data.job_types.append(jt)

            # Check basis set
            if re.match(r'\d+-\d+', kw) or 'def2' in kw.lower() or 'cc-pv' in kw.lower():
                data.basis_set = kw

            # Check SCF convergence
            for scf in scf_settings:
                if scf == kw_upper:
                    data.scf_convergence = kw

            # Check grid
            for grid in grids:
                if grid == kw_upper:
                    data.grid = kw

            # Check RI approximation
            for ri in ri_methods:
                if ri == kw_upper:
                    data.ri_approximation = kw

            # Check dispersion
            for disp in dispersion:
                if disp == kw_upper:
                    data.dispersion = kw

            # Check solvent model
            if 'CPCM(' in kw.upper():
                match = re.search(r'CPCM\(([^)]+)\)', kw, re.IGNORECASE)
                if match:
                    data.solvent = match.group(1)
            elif 'SMD' in kw.upper():
                data.solvent = 'SMD'

    # Parse % blocks
    block_pattern = r'%\s*(\w+)\s*(.*?)(?=\n%|\n\*|$)'
    for match in re.finditer(block_pattern, content, re.DOTALL | re.IGNORECASE):
        block_name = match.group(1).lower()
        block_content = match.group(2).strip()
        data.blocks[block_name] = block_content

        # Extract specific values
        if block_name == 'maxcore':
            try:
                data.maxcore = int(block_content)
            except ValueError:
                pass

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

        if data.solvent:
            print(f"Solvent: {data.solvent}")
        if data.dispersion:
            print(f"Dispersion: {data.dispersion}")
        if data.scf_convergence:
            print(f"SCF: {data.scf_convergence}")
        if data.grid:
            print(f"Grid: {data.grid}")
        if data.ri_approximation:
            print(f"RI: {data.ri_approximation}")
        if data.maxcore:
            print(f"MaxCore: {data.maxcore} MB")
        if data.nprocs > 1:
            print(f"Procs: {data.nprocs}")
