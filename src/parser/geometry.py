"""Geometry parser for ORCA output."""

import os
import re
from typing import Optional, List, Tuple
import pandas as pd

from ..core.base_parser import BaseParser
from ..core.data_models import GeometryData
from . import regex_patterns as rx


class GeometryParser(BaseParser):
    """Parse geometry data from ORCA output."""
    
    def parse(self) -> GeometryData:
        """Parse geometry data."""
        self.logger.debug("Parsing geometry...")
        
        data = GeometryData()
        
        # Parse geometry info (filename, charge, multiplicity)
        data.filename, data.charge, data.multiplicity = self._parse_geometry_info()
        
        # Parse cartesian coordinates
        data.cart_coords = self._parse_cartesian()
        if data.cart_coords is not None:
            self._log_found("atoms", len(data.cart_coords))
        
        # Generate SMILES if possible
        if data.cart_coords is not None:
            data.smiles = self._coords_to_smiles(data.cart_coords)
            if data.smiles:
                self._log_found(f"SMILES: {data.smiles}")
        
        # Parse internal coordinates
        data.internal_coords = self._parse_internal()
        
        return data
    
    def _parse_geometry_info(self) -> Tuple[Optional[str], Optional[int], Optional[int]]:
        """Parse geometry file name, charge, and multiplicity."""
        filename = None
        charge = None
        multiplicity = None
        
        # Try to find from "The coordinates will be read from file" line
        match = rx.GEOMETRY_FILE.search(self.text)
        if match:
            geometry_file = match.group(1)
            filename = os.path.splitext(geometry_file)[0]
        
        # Try to find from input section
        input_match = rx.INPUT_SECTION.search(self.text)
        if input_match:
            input_section = input_match.group(2)
            xyzfile_match = rx.XYZFILE_INPUT.search(input_section)
            if xyzfile_match:
                charge = int(xyzfile_match.group(1))
                multiplicity = int(xyzfile_match.group(2))
                geometry_file = xyzfile_match.group(3)
                filename = os.path.splitext(geometry_file)[0]
        
        return filename, charge, multiplicity
    
    def _parse_cartesian(self) -> Optional[pd.DataFrame]:
        """Parse cartesian coordinates."""
        blocks = rx.CARTESIAN_BLOCK.findall(self.text)
        if not blocks:
            self._log_not_found("cartesian coordinates")
            return None
        
        # Use last block (final geometry)
        last_block = blocks[-1].strip().splitlines()
        coords = []
        
        for line in last_block:
            parts = line.split()
            if len(parts) < 4:
                continue
            atom = parts[0]
            try:
                x = float(parts[1])
                y = float(parts[2])
                z = float(parts[3])
                coords.append({"atom": atom, "x": x, "y": y, "z": z})
            except ValueError:
                continue
        
        if not coords:
            return None
        
        return pd.DataFrame(coords)
    
    def _parse_internal(self) -> Optional[pd.DataFrame]:
        """Parse internal coordinates."""
        blocks = rx.INTERNAL_BLOCK.findall(self.text)
        if not blocks:
            return None
        
        last_block = blocks[-1].strip().splitlines()
        coords = []
        
        for line in last_block:
            parts = line.split()
            if len(parts) < 7:
                continue
            try:
                atom = parts[0]
                bond = float(parts[4])
                angle = float(parts[5])
                dihedral = float(parts[6])
                coords.append({
                    "atom": atom,
                    "bond": bond,
                    "angle": angle,
                    "dihedral": dihedral
                })
            except (ValueError, IndexError):
                continue
        
        if not coords:
            return None
        
        return pd.DataFrame(coords)
    
    def _coords_to_smiles(self, coords_df: pd.DataFrame) -> Optional[str]:
        """Convert coordinates to SMILES using OpenBabel."""
        try:
            from openbabel import pybel
            import tempfile
            
            # Write to temp XYZ file
            with tempfile.NamedTemporaryFile("w", suffix=".xyz", delete=False) as tmp:
                tmp.write(f"{len(coords_df)}\n")
                tmp.write("Generated from dataframe\n")
                for _, row in coords_df.iterrows():
                    tmp.write(f"{row['atom']} {row['x']:.6f} {row['y']:.6f} {row['z']:.6f}\n")
                tmp_path = tmp.name
            
            try:
                mol = next(pybel.readfile("xyz", tmp_path))
                smi = mol.write("smi").strip()
                return smi.split()[0] if smi else None
            finally:
                try:
                    os.unlink(tmp_path)
                except:
                    pass
                    
        except ImportError:
            self.logger.debug("OpenBabel not available, skipping SMILES generation")
            return None
        except Exception as e:
            self.logger.debug(f"SMILES generation failed: {e}")
            return None
