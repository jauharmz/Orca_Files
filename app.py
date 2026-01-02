"""ORCA Visualization Platform - Streamlit App"""

import streamlit as st
from pathlib import Path
import pandas as pd

# Page config
st.set_page_config(
    page_title="ORCA Visualization",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Initialize session state
if "df" not in st.session_state:
    st.session_state.df = None
if "parsed_files" not in st.session_state:
    st.session_state.parsed_files = []


def main():
    st.title("🔬 ORCA Quantum Chemistry Visualization")
    st.markdown("Parse and visualize ORCA output files interactively.")
    
    # Sidebar
    with st.sidebar:
        st.header("📁 Data Loading")
        
        # File upload
        uploaded_files = st.file_uploader(
            "Upload .out files",
            type=["out"],
            accept_multiple_files=True
        )
        
        if uploaded_files:
            if st.button("Parse Files"):
                parse_files(uploaded_files)
        
        # HuggingFace download
        st.divider()
        st.subheader("🤗 HuggingFace")
        if st.button("Download Sample Data"):
            download_hf_data()
        
        # Show parsed count
        if st.session_state.df is not None:
            st.success(f"✅ {len(st.session_state.df)} molecules loaded")
    
    # Main content
    if st.session_state.df is not None:
        show_data_view()
    else:
        show_welcome()


def show_welcome():
    """Show welcome message."""
    st.info("👆 Upload ORCA .out files or download sample data to get started.")
    
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("""
        ### Features
        - 📊 Energy diagrams & comparisons
        - 🧬 3D molecular viewer
        - 📈 IR, Raman, UV-Vis spectra
        - 🌳 Hierarchy & pathway detection
        - 📤 Export to JSON, CSV, HTML
        """)
    with col2:
        st.markdown("""
        ### Quick Start
        1. Upload `.out` files
        2. Click "Parse Files"
        3. Explore visualizations
        4. Export results
        """)


def parse_files(uploaded_files):
    """Parse uploaded files."""
    from src.parser.factory import ParserFactory
    
    factory = ParserFactory()
    results = []
    
    progress = st.progress(0)
    for i, file in enumerate(uploaded_files):
        try:
            content = file.read().decode('utf-8', errors='ignore')
            result = factory.parse_text(content, file.name)
            
            row = {
                "filename": file.name,
                "molecule_id": result.geometry.filename or file.name.replace('.out', ''),
                "smiles": result.geometry.smiles,
                "gibbs_Eh": result.energy.gibbs_Eh,
                "single_point_Eh": result.energy.single_point_Eh,
                "homo_energy": result.orbitals.homo_energy,
                "lumo_energy": result.orbitals.lumo_energy,
                "cart_coords": result.geometry.cart_coords,
                "orbitals": result.orbitals.orbitals,
                "ir": result.spectra.ir,
            }
            results.append(row)
        except Exception as e:
            st.warning(f"Failed to parse {file.name}: {e}")
        
        progress.progress((i + 1) / len(uploaded_files))
    
    if results:
        st.session_state.df = pd.DataFrame(results)
        st.success(f"✅ Parsed {len(results)} files")


def download_hf_data():
    """Download sample data from HuggingFace."""
    try:
        from huggingface_hub import snapshot_download
        
        with st.spinner("Downloading from HuggingFace..."):
            snapshot_download(
                repo_id="JauharMz/Orca",
                repo_type="dataset",
                local_dir="./data",
                local_dir_use_symlinks=False
            )
        st.success("✅ Downloaded to ./data")
        
        # Parse downloaded files
        files = list(Path("./data").rglob("*.out"))
        st.info(f"Found {len(files)} files. Click 'Parse Folder' to process.")
        
    except Exception as e:
        st.error(f"Download failed: {e}")


def show_data_view():
    """Show data view with tabs."""
    df = st.session_state.df
    
    # Molecule selector
    mol_ids = df["molecule_id"].dropna().unique().tolist()
    selected = st.selectbox("Select Molecule", mol_ids)
    
    # Tabs
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "📊 Energy", "🔮 Orbitals", "📈 Spectra", "📋 Data", "📤 Export"
    ])
    
    with tab1:
        show_energy_tab(df)
    
    with tab2:
        show_orbital_tab(df, selected)
    
    with tab3:
        show_spectra_tab(df, selected)
    
    with tab4:
        show_data_tab(df)
    
    with tab5:
        show_export_tab(df)


def show_energy_tab(df):
    """Energy comparison view."""
    st.subheader("Energy Comparison")
    
    if "gibbs_Eh" not in df.columns or df["gibbs_Eh"].isna().all():
        st.warning("No energy data available")
        return
    
    from src.viz.energy_diagram import EnergyDiagramVisualizer
    
    viz = EnergyDiagramVisualizer(df)
    fig = viz.create_figure()
    st.plotly_chart(fig, use_container_width=True)


def show_orbital_tab(df, selected):
    """Orbital view."""
    st.subheader(f"Orbitals - {selected}")
    
    row = df[df["molecule_id"] == selected]
    if row.empty or "orbitals" not in row.columns:
        st.warning("No orbital data")
        return
    
    orbitals = row.iloc[0]["orbitals"]
    if orbitals is None or orbitals.empty:
        st.warning("No orbital data")
        return
    
    from src.viz.orbital_plot import OrbitalPlotVisualizer
    
    viz = OrbitalPlotVisualizer(orbitals)
    fig = viz.create_figure()
    st.plotly_chart(fig, use_container_width=True)


def show_spectra_tab(df, selected):
    """Spectra view."""
    st.subheader(f"Spectra - {selected}")
    
    row = df[df["molecule_id"] == selected]
    if row.empty:
        return
    
    row_data = row.iloc[0]
    
    if "ir" in row_data and row_data["ir"] is not None:
        from src.viz.spectra_plot import SpectraVisualizer
        viz = SpectraVisualizer(row_data["ir"])
        fig = viz.create_ir_figure()
        st.plotly_chart(fig, use_container_width=True)
    else:
        st.info("No IR spectrum data")


def show_data_tab(df):
    """Raw data view."""
    st.subheader("Parsed Data")
    
    # Show scalar columns
    scalar_cols = ["molecule_id", "smiles", "gibbs_Eh", "single_point_Eh", 
                   "homo_energy", "lumo_energy"]
    available = [c for c in scalar_cols if c in df.columns]
    
    st.dataframe(df[available], use_container_width=True)


def show_export_tab(df):
    """Export view."""
    st.subheader("Export Data")
    
    col1, col2 = st.columns(2)
    
    with col1:
        if st.button("📥 Export JSON"):
            from src.export.data_exporter import DataExporter
            exporter = DataExporter(df)
            path = exporter.to_json("export.json")
            st.success(f"Exported: {path}")
    
    with col2:
        if st.button("📥 Export CSV"):
            from src.export.data_exporter import DataExporter
            exporter = DataExporter(df)
            path = exporter.to_csv("export.csv")
            st.success(f"Exported: {path}")
    
    if st.button("📄 Export HTML Report"):
        from src.export.html_exporter import HTMLExporter
        exporter = HTMLExporter(df)
        path = exporter.export("report.html")
        st.success(f"Exported: {path}")


if __name__ == "__main__":
    main()
