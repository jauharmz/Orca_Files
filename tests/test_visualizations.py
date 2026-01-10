"""
Visualization Test - Test all visualizations with parsed data

Creates HTML files with interactive plots.

Run: python tests/test_visualizations.py
"""

import sys
from pathlib import Path
import logging
import pandas as pd

# Add root to path
ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(name)-25s | %(message)s'
)
logger = logging.getLogger("test_visualizations")


def main():
    logger.info("=" * 70)
    logger.info("VISUALIZATION TEST")
    logger.info("=" * 70)
    
    # 1. Load data
    logger.info("\n[STEP 1] Loading data...")
    from huggingface_hub import snapshot_download
    
    data_dir = Path("./test_data_hf")
    if not data_dir.exists():
        snapshot_download(
            repo_id="JauharMz/Orca",
            repo_type="dataset",
            local_dir=str(data_dir),
            local_dir_use_symlinks=False
        )
    
    # Parse data
    from src.parser.batch import BatchParser
    pattern = str(data_dir / "**/*.out")
    batch = BatchParser(pattern)
    df = batch.parse_all(verbose=False)
    
    logger.info(f"  Loaded {len(df)} molecules")
    
    # 2. Energy diagram
    logger.info("\n[STEP 2] Creating energy diagram...")
    from src.viz.energy_diagram import EnergyDiagramVisualizer
    
    energy_viz = EnergyDiagramVisualizer(df, id_column="molecule_id")
    
    # Try gibbs first, then single point
    if df["gibbs_Eh"].notna().any():
        fig_energy = energy_viz.create_figure(energy_col="gibbs_Eh")
    else:
        fig_energy = energy_viz.create_figure(energy_col="single_point_Eh")
    
    fig_energy.write_html("viz_energy.html")
    logger.info("  Saved: viz_energy.html")
    
    # 3. Pathway diagram (p1x -> p2x -> p3x -> ...)
    logger.info("\n[STEP 3] Creating pathway diagram...")
    pathway = ["p1xs0", "p2xs0", "p3xs0", "p4xs0", "p5xs0", "p6xs0"]
    available = df[df["molecule_id"].isin(pathway)]["molecule_id"].tolist()
    
    if len(available) >= 2:
        fig_pathway = energy_viz.create_pathway_figure(available, energy_col="gibbs_Eh")
        fig_pathway.write_html("viz_pathway.html")
        logger.info(f"  Saved: viz_pathway.html ({len(available)} steps)")
    else:
        logger.warning("  Not enough molecules for pathway")
    
    # 4. Orbital plot (first molecule with orbital data)
    logger.info("\n[STEP 4] Creating orbital plot...")
    from src.viz.orbital_plot import OrbitalPlotVisualizer
    
    # Find a molecule with orbital data
    orbital_data = None
    for _, row in df.iterrows():
        orbitals = row.get("orbitals")
        if isinstance(orbitals, pd.DataFrame) and not orbitals.empty:
            orbital_data = orbitals
            logger.info(f"  Using orbitals from: {row['molecule_id']}")
            break
    
    if orbital_data is not None:
        orb_viz = OrbitalPlotVisualizer(orbital_data)
        fig_orb = orb_viz.create_figure(n_orbitals=5)
        fig_orb.write_html("viz_orbitals.html")
        logger.info("  Saved: viz_orbitals.html")
    else:
        logger.warning("  No orbital data found")
    
    # 5. IR spectrum (first molecule with IR data)
    logger.info("\n[STEP 5] Creating IR spectrum...")
    from src.viz.spectra_plot import SpectraVisualizer
    
    ir_data = None
    for _, row in df.iterrows():
        ir = row.get("ir")
        if isinstance(ir, pd.DataFrame) and not ir.empty:
            ir_data = ir
            logger.info(f"  Using IR from: {row['molecule_id']}")
            break
    
    if ir_data is not None:
        ir_viz = SpectraVisualizer(ir_data)
        fig_ir = ir_viz.create_ir_figure()
        fig_ir.write_html("viz_ir.html")
        logger.info("  Saved: viz_ir.html")
    else:
        logger.warning("  No IR data found")
    
    # 6. UV-Vis spectrum (from TDDFT data)
    logger.info("\n[STEP 6] Creating UV-Vis spectrum...")
    
    uv_data = None
    for _, row in df.iterrows():
        tddft = row.get("tddft_states")
        if isinstance(tddft, pd.DataFrame) and not tddft.empty:
            # Create UV-Vis from TDDFT states
            if "energy_ev" in tddft.columns:
                # Convert energy to wavelength (nm)
                tddft = tddft.copy()
                tddft["wavelength_nm"] = 1239.84 / tddft["energy_ev"]
                tddft["fosc"] = tddft.get("weight", 1.0).abs()
                uv_data = tddft[["wavelength_nm", "fosc"]].drop_duplicates()
                logger.info(f"  Using TDDFT from: {row['molecule_id']}")
                break
    
    if uv_data is not None:
        uv_viz = SpectraVisualizer(uv_data)
        fig_uv = uv_viz.create_uv_vis_figure()
        fig_uv.write_html("viz_uvvis.html")
        logger.info("  Saved: viz_uvvis.html")
    else:
        logger.warning("  No UV-Vis data found")
    
    # 7. Summary
    logger.info("\n" + "=" * 70)
    logger.info("VISUALIZATION TEST COMPLETE!")
    logger.info("Created files:")
    logger.info("  - viz_energy.html (energy bar chart)")
    logger.info("  - viz_pathway.html (energy pathway)")
    logger.info("  - viz_orbitals.html (orbital levels)")
    logger.info("  - viz_ir.html (IR spectrum)")
    logger.info("  - viz_uvvis.html (UV-Vis spectrum)")
    logger.info("=" * 70)
    
    return df


if __name__ == "__main__":
    df = main()
