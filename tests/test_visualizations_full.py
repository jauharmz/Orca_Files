"""
Comprehensive Visualization Test

Creates complete set of visualizations including:
- Energy diagrams
- Pathway diagrams
- Orbital plots
- IR and Raman spectra
- UV-Vis spectra
- IR-Raman comparison (dual axis)
- Multi-molecule energy comparison
- 3D molecular structure (py3Dmol)

Run: python tests/test_visualizations_full.py
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
logger = logging.getLogger("test_viz_full")

# Output directory
OUTPUT_DIR = Path("./viz_output")
OUTPUT_DIR.mkdir(exist_ok=True)


def load_data():
    """Load and parse data from HuggingFace."""
    logger.info("=" * 70)
    logger.info("LOADING DATA")
    logger.info("=" * 70)
    
    from huggingface_hub import snapshot_download
    
    data_dir = Path("./test_data_hf")
    if not data_dir.exists():
        snapshot_download(
            repo_id="JauharMz/Orca",
            repo_type="dataset",
            local_dir=str(data_dir),
            local_dir_use_symlinks=False
        )
    
    from src.parser.batch import BatchParser
    pattern = str(data_dir / "**/*.out")
    batch = BatchParser(pattern)
    df = batch.parse_all(verbose=False)
    
    logger.info(f"  Loaded {len(df)} molecules")
    return df


def create_energy_diagram(df):
    """Create energy level diagram."""
    logger.info("\n[1] ENERGY DIAGRAM")
    
    from src.viz.energy_diagram import EnergyDiagramVisualizer
    
    viz = EnergyDiagramVisualizer(df, id_column="molecule_id")
    
    if df["gibbs_Eh"].notna().any():
        fig = viz.create_figure(energy_col="gibbs_Eh")
        fig.update_layout(title="Gibbs Free Energy (Eh)")
    else:
        fig = viz.create_figure(energy_col="single_point_Eh")
        fig.update_layout(title="Single Point Energy (Eh)")
    
    path = OUTPUT_DIR / "energy_diagram.html"
    fig.write_html(str(path))
    logger.info(f"  Saved: {path}")
    return fig


def create_pathway_diagram(df):
    """Create energy pathway diagram."""
    logger.info("\n[2] PATHWAY DIAGRAM")
    
    from src.viz.energy_diagram import EnergyDiagramVisualizer
    
    # Try to find pathway molecules
    pathway_patterns = [
        ["p1xs0", "p2xs0", "p3xs0", "p4xs0", "p5xs0", "p6xs0"],
        ["p1x", "p2x", "p3x", "p4x", "p5x", "p6x"],
    ]
    
    available = []
    for pattern in pathway_patterns:
        available = df[df["molecule_id"].isin(pattern)]["molecule_id"].tolist()
        if len(available) >= 2:
            break
    
    if len(available) >= 2:
        viz = EnergyDiagramVisualizer(df, id_column="molecule_id")
        fig = viz.create_pathway_figure(available, energy_col="gibbs_Eh")
        fig.update_layout(title=f"Energy Pathway ({len(available)} steps)")
        
        path = OUTPUT_DIR / "pathway_diagram.html"
        fig.write_html(str(path))
        logger.info(f"  Saved: {path}")
        return fig
    else:
        logger.warning("  Not enough molecules for pathway")
        return None


def create_orbital_plot(df):
    """Create orbital energy level plot."""
    logger.info("\n[3] ORBITAL PLOT")
    
    from src.viz.orbital_plot import OrbitalPlotVisualizer
    
    for _, row in df.iterrows():
        orbitals = row.get("orbitals")
        if isinstance(orbitals, pd.DataFrame) and not orbitals.empty:
            logger.info(f"  Using: {row['molecule_id']}")
            
            viz = OrbitalPlotVisualizer(orbitals)
            fig = viz.create_figure(n_orbitals=10)
            fig.update_layout(title=f"Orbital Energy Levels - {row['molecule_id']}")
            
            path = OUTPUT_DIR / "orbitals.html"
            fig.write_html(str(path))
            logger.info(f"  Saved: {path}")
            return fig
    
    logger.warning("  No orbital data found")
    return None


def create_ir_spectrum(df):
    """Create IR spectrum."""
    logger.info("\n[4] IR SPECTRUM")
    
    from src.viz.spectra_plot import SpectraVisualizer
    
    for _, row in df.iterrows():
        ir = row.get("ir")
        if isinstance(ir, pd.DataFrame) and not ir.empty:
            logger.info(f"  Using: {row['molecule_id']} ({len(ir)} peaks)")
            
            viz = SpectraVisualizer(ir)
            fig = viz.create_ir_figure()
            fig.update_layout(title=f"IR Spectrum - {row['molecule_id']}")
            
            path = OUTPUT_DIR / "ir_spectrum.html"
            fig.write_html(str(path))
            logger.info(f"  Saved: {path}")
            return fig, ir, row['molecule_id']
    
    logger.warning("  No IR data found")
    return None, None, None


def create_raman_spectrum(df):
    """Create Raman spectrum."""
    logger.info("\n[5] RAMAN SPECTRUM")
    
    from src.viz.spectra_plot import SpectraVisualizer
    
    for _, row in df.iterrows():
        raman = row.get("raman")
        if isinstance(raman, pd.DataFrame) and not raman.empty:
            logger.info(f"  Using: {row['molecule_id']} ({len(raman)} peaks)")
            
            viz = SpectraVisualizer(raman)
            fig = viz.create_raman_figure()
            fig.update_layout(title=f"Raman Spectrum - {row['molecule_id']}")
            
            path = OUTPUT_DIR / "raman_spectrum.html"
            fig.write_html(str(path))
            logger.info(f"  Saved: {path}")
            return fig, raman, row['molecule_id']
    
    logger.warning("  No Raman data found")
    return None, None, None


def create_ir_raman_comparison(ir_data, raman_data, mol_id):
    """Create IR-Raman comparison with dual-axis."""
    logger.info("\n[6] IR-RAMAN COMPARISON")
    
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    
    if ir_data is None or raman_data is None:
        logger.warning("  Need both IR and Raman data")
        return None
    
    # Create dual-axis plot
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    
    # IR trace
    fig.add_trace(
        go.Scatter(
            x=ir_data["freq_cm-1"], 
            y=ir_data["intensity_km/mol"],
            mode="lines",
            name="IR",
            line=dict(color="blue", width=1.5),
            fill="tozeroy",
            fillcolor="rgba(0,0,255,0.1)"
        ),
        secondary_y=False
    )
    
    # Raman trace
    fig.add_trace(
        go.Scatter(
            x=raman_data["freq_cm-1"],
            y=raman_data["activity"],
            mode="lines",
            name="Raman",
            line=dict(color="green", width=1.5),
            fill="tozeroy",
            fillcolor="rgba(0,255,0,0.1)"
        ),
        secondary_y=True
    )
    
    fig.update_layout(
        title=f"IR-Raman Comparison - {mol_id}",
        xaxis_title="Frequency (cm⁻¹)",
        xaxis=dict(autorange="reversed"),
        hovermode="x unified",
        legend=dict(x=0.8, y=0.95)
    )
    fig.update_yaxes(title_text="IR Intensity (km/mol)", secondary_y=False, color="blue")
    fig.update_yaxes(title_text="Raman Activity", secondary_y=True, color="green")
    
    path = OUTPUT_DIR / "ir_raman_comparison.html"
    fig.write_html(str(path))
    logger.info(f"  Saved: {path}")
    return fig


def create_uv_vis_spectrum(df):
    """Create UV-Vis spectrum from TDDFT data."""
    logger.info("\n[7] UV-VIS SPECTRUM")
    
    import plotly.graph_objects as go
    import numpy as np
    
    for _, row in df.iterrows():
        tddft = row.get("tddft_states")
        if isinstance(tddft, pd.DataFrame) and not tddft.empty:
            if "energy_ev" not in tddft.columns:
                continue
            
            logger.info(f"  Using: {row['molecule_id']} ({len(tddft)} states)")
            
            # Convert eV to nm
            tddft = tddft.copy()
            tddft["wavelength_nm"] = 1239.84 / tddft["energy_ev"]
            tddft["fosc"] = tddft.get("weight", pd.Series([0.5]*len(tddft))).abs()
            
            # Create stick spectrum
            fig = go.Figure()
            
            # Stick plot
            for _, state in tddft.iterrows():
                fig.add_trace(go.Scatter(
                    x=[state["wavelength_nm"], state["wavelength_nm"]],
                    y=[0, state["fosc"]],
                    mode="lines",
                    line=dict(color="purple", width=2),
                    showlegend=False,
                    hovertemplate=f"λ={state['wavelength_nm']:.1f}nm<br>f={state['fosc']:.3f}"
                ))
            
            # Generate broadened spectrum
            wl_range = np.linspace(200, 800, 600)
            sigma = 20  # nm
            spectrum = np.zeros_like(wl_range)
            for _, state in tddft.iterrows():
                spectrum += state["fosc"] * np.exp(-0.5 * ((wl_range - state["wavelength_nm"]) / sigma) ** 2)
            
            fig.add_trace(go.Scatter(
                x=wl_range,
                y=spectrum,
                mode="lines",
                name="Broadened",
                line=dict(color="rgba(128,0,128,0.5)", width=1),
                fill="tozeroy",
                fillcolor="rgba(128,0,128,0.1)"
            ))
            
            fig.update_layout(
                title=f"UV-Vis Spectrum - {row['molecule_id']}",
                xaxis_title="Wavelength (nm)",
                yaxis_title="Oscillator Strength",
                xaxis=dict(range=[200, 800]),
                showlegend=False
            )
            
            path = OUTPUT_DIR / "uvvis_spectrum.html"
            fig.write_html(str(path))
            logger.info(f"  Saved: {path}")
            return fig
    
    logger.warning("  No TDDFT data found")
    return None


def create_multi_molecule_comparison(df):
    """Create multi-molecule energy comparison."""
    logger.info("\n[8] MULTI-MOLECULE COMPARISON")
    
    import plotly.graph_objects as go
    
    # Get molecules with Gibbs energy
    subset = df[df["gibbs_Eh"].notna()].head(10)
    
    if len(subset) < 2:
        logger.warning("  Not enough molecules with energy data")
        return None
    
    # Convert to relative energy (kcal/mol)
    min_energy = subset["gibbs_Eh"].min()
    subset = subset.copy()
    subset["rel_energy_kcal"] = (subset["gibbs_Eh"] - min_energy) * 627.509
    
    fig = go.Figure()
    
    # Bar chart
    colors = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', 
              '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']
    
    fig.add_trace(go.Bar(
        x=subset["molecule_id"],
        y=subset["rel_energy_kcal"],
        marker_color=colors[:len(subset)],
        text=[f"{e:.1f}" for e in subset["rel_energy_kcal"]],
        textposition="outside"
    ))
    
    fig.update_layout(
        title="Relative Energy Comparison",
        xaxis_title="Molecule",
        yaxis_title="Relative Energy (kcal/mol)",
        yaxis=dict(range=[0, subset["rel_energy_kcal"].max() * 1.2])
    )
    
    path = OUTPUT_DIR / "energy_comparison.html"
    fig.write_html(str(path))
    logger.info(f"  Saved: {path}")
    return fig


def create_3d_structure(df):
    """Create 3D molecular structure using py3Dmol."""
    logger.info("\n[9] 3D MOLECULAR STRUCTURE")
    
    try:
        import py3Dmol
    except ImportError:
        logger.warning("  py3Dmol not installed. Run: pip install py3Dmol")
        return None
    
    for _, row in df.iterrows():
        coords = row.get("cart_coords")
        if isinstance(coords, pd.DataFrame) and not coords.empty:
            mol_id = row['molecule_id']
            logger.info(f"  Using: {mol_id} ({len(coords)} atoms)")
            
            # Create XYZ format
            xyz_lines = [str(len(coords)), mol_id]
            for _, atom in coords.iterrows():
                xyz_lines.append(f"{atom['atom']}  {atom['x']:.6f}  {atom['y']:.6f}  {atom['z']:.6f}")
            xyz_str = "\n".join(xyz_lines)
            
            # Create py3Dmol viewer
            viewer = py3Dmol.view(width=800, height=600)
            viewer.addModel(xyz_str, "xyz")
            viewer.setStyle({"stick": {"radius": 0.1}, "sphere": {"scale": 0.25}})
            viewer.setBackgroundColor("white")
            viewer.zoomTo()
            
            # Get HTML
            html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>3D Structure - {mol_id}</title>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <style>
        body {{ margin: 0; padding: 20px; font-family: Arial, sans-serif; }}
        h1 {{ color: #333; }}
        #viewer {{ width: 800px; height: 600px; border: 1px solid #ddd; }}
    </style>
</head>
<body>
    <h1>3D Molecular Structure - {mol_id}</h1>
    <p>Atoms: {len(coords)} | SMILES: {row.get('smiles', 'N/A')}</p>
    <div id="viewer"></div>
    <script>
        var viewer = $3Dmol.createViewer("viewer", {{backgroundColor: "white"}});
        var xyz = `{xyz_str}`;
        viewer.addModel(xyz, "xyz");
        viewer.setStyle({{}}, {{stick: {{radius: 0.15}}, sphere: {{scale: 0.3}}}});
        viewer.zoomTo();
        viewer.render();
    </script>
</body>
</html>
"""
            
            path = OUTPUT_DIR / "structure_3d.html"
            with open(path, "w") as f:
                f.write(html_content)
            logger.info(f"  Saved: {path}")
            return path
    
    logger.warning("  No coordinate data found")
    return None


def create_orbital_comparison(df):
    """Create HOMO-LUMO comparison across molecules."""
    logger.info("\n[10] HOMO-LUMO COMPARISON")
    
    import plotly.graph_objects as go
    
    subset = df[df["homo_energy"].notna() & df["lumo_energy"].notna()].head(10)
    
    if len(subset) < 2:
        logger.warning("  Not enough molecules with orbital data")
        return None
    
    fig = go.Figure()
    
    # HOMO bars
    fig.add_trace(go.Bar(
        name="HOMO",
        x=subset["molecule_id"],
        y=subset["homo_energy"],
        marker_color="blue"
    ))
    
    # LUMO bars
    fig.add_trace(go.Bar(
        name="LUMO",
        x=subset["molecule_id"],
        y=subset["lumo_energy"],
        marker_color="red"
    ))
    
    # Gap annotations
    for _, row in subset.iterrows():
        gap = row["lumo_energy"] - row["homo_energy"]
        fig.add_annotation(
            x=row["molecule_id"],
            y=(row["homo_energy"] + row["lumo_energy"]) / 2,
            text=f"Δ={gap:.2f}eV",
            showarrow=False,
            font=dict(size=10)
        )
    
    fig.update_layout(
        title="HOMO-LUMO Comparison",
        xaxis_title="Molecule",
        yaxis_title="Energy (eV)",
        barmode="group",
        legend=dict(x=0.8, y=0.95)
    )
    
    path = OUTPUT_DIR / "homo_lumo_comparison.html"
    fig.write_html(str(path))
    logger.info(f"  Saved: {path}")
    return fig


def main():
    """Run all visualization tests."""
    logger.info("=" * 70)
    logger.info("COMPREHENSIVE VISUALIZATION TEST")
    logger.info("=" * 70)
    
    # Load data
    df = load_data()
    
    # Create visualizations
    create_energy_diagram(df)
    create_pathway_diagram(df)
    create_orbital_plot(df)
    ir_fig, ir_data, ir_mol = create_ir_spectrum(df)
    raman_fig, raman_data, raman_mol = create_raman_spectrum(df)
    
    # IR-Raman comparison (if both exist)
    if ir_data is not None and raman_data is not None:
        create_ir_raman_comparison(ir_data, raman_data, ir_mol)
    
    create_uv_vis_spectrum(df)
    create_multi_molecule_comparison(df)
    create_3d_structure(df)
    create_orbital_comparison(df)
    
    # Summary
    logger.info("\n" + "=" * 70)
    logger.info("VISUALIZATION TEST COMPLETE!")
    logger.info(f"Output directory: {OUTPUT_DIR.absolute()}")
    logger.info("Created files:")
    for f in sorted(OUTPUT_DIR.glob("*.html")):
        logger.info(f"  - {f.name}")
    logger.info("=" * 70)
    
    return df


if __name__ == "__main__":
    df = main()
