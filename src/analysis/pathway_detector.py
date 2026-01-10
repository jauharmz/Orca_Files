"""Pathway detector for degradation pathways."""

import re
from typing import Dict, List, Optional, Tuple, Any
from collections import defaultdict
import pandas as pd

from ..core.data_models import Pathway, Hierarchy
from ..logger import get_logger


class PathwayDetector:
    """
    Detect degradation pathways from molecule naming patterns.
    
    Supports:
        - Automatic edge detection from naming patterns
        - Reaction rules (add/remove species)
        - Step corrections for branch initialization
        - Color schemes (by_variant, by_destination)
    """
    
    def __init__(
        self,
        df: pd.DataFrame,
        hierarchy: Optional[Hierarchy] = None,
        id_column: str = "molecule_id"
    ):
        """
        Initialize detector.
        
        Args:
            df: DataFrame with molecule data
            hierarchy: Optional pre-computed hierarchy
            id_column: Column containing molecule IDs
        """
        self.df = df
        self.hierarchy = hierarchy
        self.id_column = id_column
        self.logger = get_logger("PathwayDetector")
        
        # Configuration
        self.reaction_rules: Dict[Tuple[str, str], Dict] = {}
        self.step_corrections: Dict[Tuple[str, str], Dict] = {}
        self.color_scheme = "by_variant"
        self.colors: Dict[str, List[List[str]]] = {}
    
    def set_reaction_rules(self, rules: Dict[Tuple[str, str], Dict]):
        """
        Set reaction rules for pathway transitions.
        
        Args:
            rules: Dict mapping (from_pattern, to_pattern) to reaction
                   e.g., {("p1", "p2"): {"add": {"OH": 4}, "remove": {"H2O": 3}}}
        """
        self.reaction_rules = rules
        self.logger.info(f"Set {len(rules)} reaction rules")
    
    def set_step_corrections(self, corrections: Dict[Tuple[str, str], Dict]):
        """
        Set step corrections for branch initialization.
        
        Args:
            corrections: Dict mapping (from_id, to_id) to correction
                        e.g., {("p1x", "p1a"): {"add": {"OH": 2}}}
        """
        self.step_corrections = corrections
        self.logger.info(f"Set {len(corrections)} step corrections")
    
    def set_color_scheme(self, scheme: str):
        """Set color scheme: 'by_variant' or 'by_destination'."""
        self.color_scheme = scheme
    
    def set_colors(self, colors: Dict[str, List[List[str]]]):
        """
        Set custom color scheme.
        
        Args:
            colors: Dict mapping color to list of pathway segments
                   e.g., {"blue": [["p1x", "p2x", "p3x"]]}
        """
        self.colors = colors
    
    def detect(self) -> List[Pathway]:
        """
        Detect all pathways.
        
        Returns:
            List of Pathway objects
        """
        if self.id_column not in self.df.columns:
            self.logger.warning(f"Column {self.id_column} not found")
            return []
        
        ids = self.df[self.id_column].dropna().unique().tolist()
        self.logger.info(f"Detecting pathways from {len(ids)} molecules")
        
        # Detect edges
        edges = self._detect_edges(ids)
        self.logger.info(f"Found {len(edges)} edges")
        
        # Build pathways
        pathways = self._build_pathways(edges)
        
        # Apply reaction rules
        for pathway in pathways:
            pathway.reactions = self._get_reactions(pathway.edges)
        
        # Apply colors
        self._apply_colors(pathways)
        
        return pathways
    
    def _detect_edges(self, ids: List[str]) -> List[Tuple[str, str]]:
        """Detect edges from molecule naming patterns."""
        edges = []
        id_set = set(ids)
        
        for mol_id in ids:
            # Extract base pattern and step number
            match = re.match(r'^(.+?)(\d+)([a-z]?)$', mol_id, re.IGNORECASE)
            if not match:
                continue
            
            prefix = match.group(1)
            step = int(match.group(2))
            variant = match.group(3).lower() if match.group(3) else ""
            
            # Look for next step
            next_id = f"{prefix}{step + 1}{variant}"
            if next_id in id_set:
                edges.append((mol_id, next_id))
            
            # Look for branches (same prefix + different variant)
            for other_id in ids:
                other_match = re.match(r'^(.+?)(\d+)([a-z]?)$', other_id, re.IGNORECASE)
                if not other_match:
                    continue
                
                other_prefix = other_match.group(1)
                other_step = int(other_match.group(2))
                other_variant = other_match.group(3).lower() if other_match.group(3) else ""
                
                # Branch from main to variant
                if prefix == other_prefix and step == other_step and variant != other_variant:
                    if not variant and other_variant:  # main → variant
                        edges.append((mol_id, other_id))
        
        return list(set(edges))
    
    def _build_pathways(self, edges: List[Tuple[str, str]]) -> List[Pathway]:
        """Build pathway objects from edges."""
        # Group edges by variant
        variant_edges = defaultdict(list)
        
        for from_id, to_id in edges:
            # Extract variant
            match = re.search(r'([a-z])$', from_id, re.IGNORECASE)
            variant = match.group(1).lower() if match else "main"
            variant_edges[variant].append((from_id, to_id))
        
        pathways = []
        for variant, var_edges in variant_edges.items():
            # Get all nodes
            nodes = set()
            for from_id, to_id in var_edges:
                nodes.add(from_id)
                nodes.add(to_id)
            
            # Sort nodes
            def sort_key(n):
                match = re.search(r'(\d+)', n)
                return int(match.group(1)) if match else 0
            
            sorted_nodes = sorted(nodes, key=sort_key)
            
            pathway = Pathway(
                nodes=sorted_nodes,
                edges=var_edges,
                color=None
            )
            pathways.append(pathway)
        
        return pathways
    
    def _get_reactions(self, edges: List[Tuple[str, str]]) -> Dict[Tuple[str, str], Dict]:
        """Get reactions for pathway edges."""
        reactions = {}
        
        for from_id, to_id in edges:
            # Check direct match
            if (from_id, to_id) in self.step_corrections:
                reactions[(from_id, to_id)] = self.step_corrections[(from_id, to_id)]
                continue
            
            # Check pattern match
            for (from_pattern, to_pattern), reaction in self.reaction_rules.items():
                if from_pattern in from_id and to_pattern in to_id:
                    reactions[(from_id, to_id)] = reaction
                    break
        
        return reactions
    
    def _apply_colors(self, pathways: List[Pathway]):
        """Apply color scheme to pathways."""
        if self.colors:
            # Use custom colors
            for color, segments in self.colors.items():
                for pathway in pathways:
                    for segment in segments:
                        if all(n in pathway.nodes for n in segment):
                            pathway.color = color
                            break
        else:
            # Auto-assign colors
            colors = ["blue", "red", "green", "orange", "purple", "brown", "pink"]
            for i, pathway in enumerate(pathways):
                pathway.color = colors[i % len(colors)]
    
    def to_edges_list(self) -> List[Tuple[str, str]]:
        """Get all edges as list."""
        pathways = self.detect()
        edges = []
        for p in pathways:
            edges.extend(p.edges)
        return list(set(edges))
