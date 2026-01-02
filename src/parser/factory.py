"""
Parser Factory - Orchestrates all modular parsers

This is the main entry point for parsing ORCA output files.
"""

import os
from pathlib import Path
from typing import Optional, Dict, Any
import pandas as pd

from .geometry import GeometryParser
from .energy import EnergyParser
from .orbitals import OrbitalParser
from .spectroscopy import SpectroscopyParser
from .tddft import TDDFTParser
from ..core.data_models import ParseResult
from ..logger import get_logger


class ParserFactory:
    """Factory that orchestrates all modular parsers."""
    
    def __init__(self):
        self.logger = get_logger("ParserFactory")
    
    def parse(self, filepath: str) -> ParseResult:
        """Parse an ORCA output file."""
        self.logger.info(f"Parsing file: {filepath}")
        
        try:
            with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
                text = f.read()
            self.logger.debug(f"  File size: {len(text)} bytes")
            return self._parse_text(text, filepath)
        except Exception as e:
            self.logger.error(f"Parse failed: {e}")
            return self._empty_result(filepath)
    
    def parse_text(self, text: str, filename: str = "unknown") -> ParseResult:
        """Parse from text content."""
        self.logger.info(f"Parsing text: {filename} ({len(text)} bytes)")
        return self._parse_text(text, filename)
    
    def _parse_text(self, text: str, filename: str) -> ParseResult:
        """Internal parsing logic."""
        if not text:
            self.logger.warning("Empty text, returning empty result")
            return self._empty_result(filename)
        
        self.logger.debug("=" * 50)
        self.logger.debug(f"PARSING: {filename}")
        self.logger.debug("=" * 50)
        
        # Run all parsers
        self.logger.debug("--- Geometry Parser ---")
        geo_parser = GeometryParser(text)
        geometry = geo_parser.parse()
        internal_coords = geo_parser.parse_internal()
        
        self.logger.debug("--- Energy Parser ---")
        energy_parser = EnergyParser(text)
        energy = energy_parser.parse()
        calc_info = energy_parser.detect_calc_type()
        
        self.logger.debug("--- Orbital Parser ---")
        orbital_parser = OrbitalParser(text)
        orbitals = orbital_parser.parse()
        
        self.logger.debug("--- Spectroscopy Parser ---")
        spec_parser = SpectroscopyParser(text)
        spectra = spec_parser.parse()
        mulliken = spec_parser.parse_mulliken()
        
        self.logger.debug("--- TD-DFT Parser ---")
        tddft_parser = TDDFTParser(text)
        tddft = tddft_parser.parse()
        
        # Summary
        self.logger.debug("=" * 50)
        self.logger.info(f"Parse complete: {filename}")
        self._log_summary(geometry, energy, orbitals, spectra, tddft, calc_info)
        
        return ParseResult(
            filename=filename,
            geometry=geometry,
            energy=energy,
            orbitals=orbitals,
            spectra=spectra,
            tddft=tddft,
            mulliken=mulliken,
            internal_coords=internal_coords,
            is_optimization=calc_info.get("is_optimization", False),
            has_tddft=calc_info.get("has_tddft", False),
            optimized_state=calc_info.get("optimized_state", "S0"),
            esd_type=calc_info.get("esd_type"),
            calc_class=calc_info.get("calc_class", "single_point")
        )
    
    def _log_summary(self, geometry, energy, orbitals, spectra, tddft, calc_info):
        """Log parsing summary."""
        summary = []
        
        if geometry.filename:
            summary.append(f"mol={geometry.filename}")
        if geometry.smiles:
            summary.append(f"smiles=yes")
        if geometry.cart_coords is not None:
            summary.append(f"atoms={len(geometry.cart_coords)}")
        if energy.gibbs_Eh:
            summary.append(f"gibbs={energy.gibbs_Eh:.4f}")
        if energy.single_point_Eh:
            summary.append(f"sp={energy.single_point_Eh:.4f}")
        if orbitals.homo_energy:
            summary.append(f"homo={orbitals.homo_energy:.2f}eV")
        if orbitals.lumo_energy:
            summary.append(f"lumo={orbitals.lumo_energy:.2f}eV")
        if spectra.ir is not None and not spectra.ir.empty:
            summary.append(f"ir={len(spectra.ir)}peaks")
        if tddft.states is not None and not tddft.states.empty:
            summary.append(f"tddft={len(tddft.states)}trans")
        
        summary.append(f"class={calc_info.get('calc_class')}")
        
        self.logger.info(f"  Summary: {', '.join(summary)}")
    
    def _empty_result(self, filepath: str) -> ParseResult:
        """Return empty result."""
        from ..core.data_models import GeometryData, EnergyData, OrbitalData, SpectraData, TDDFTData, MullikenData, InternalCoordsData
        
        return ParseResult(
            filename=filepath,
            geometry=GeometryData(),
            energy=EnergyData(),
            orbitals=OrbitalData(),
            spectra=SpectraData(),
            tddft=TDDFTData(),
            mulliken=MullikenData(),
            internal_coords=InternalCoordsData()
        )
