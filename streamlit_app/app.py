"""
ORCA Visualization Platform v2 - Main Entry

Run: streamlit run streamlit_app/app.py
"""

import streamlit as st
import sys
from pathlib import Path

# Add parent to path for src imports
ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

# Page config
st.set_page_config(
    page_title="ORCA Visualization",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Import utilities
from utils.session_state import init_session_state, get_df, set_df

# Initialize session state
init_session_state()


def main():
    """Main app entry."""
    
    # Sidebar - Data Loading
    with st.sidebar:
        st.title("🔬 ORCA Platform")
        st.divider()
        
        # Data loading section
        st.header("📁 Data Loading")
        
        data_source = st.radio(
            "Source",
            ["Upload Files", "HuggingFace", "Local Folder"],
            horizontal=True
        )
        
        if data_source == "Upload Files":
            uploaded = st.file_uploader(
                "Upload .out files",
                type=["out"],
                accept_multiple_files=True
            )
            if uploaded and st.button("🚀 Parse Files", type="primary"):
                parse_uploaded_files(uploaded)
        
        elif data_source == "HuggingFace":
            # Check if data already exists
            data_path = Path("./test_data_hf")
            if data_path.exists() and list(data_path.rglob("*.out")):
                st.info(f"📁 Sample data exists ({len(list(data_path.rglob('*.out')))} files)")
                col1, col2 = st.columns(2)
                with col1:
                    if st.button("🚀 Parse Existing", type="primary"):
                        parse_folder("./test_data_hf")
                with col2:
                    if st.button("📥 Re-download"):
                        download_hf_data()
            else:
                if st.button("📥 Download & Parse", type="primary"):
                    download_hf_data()
        
        elif data_source == "Local Folder":
            folder = st.text_input("Folder path", "./test_data_hf")
            if st.button("📂 Parse Folder", type="primary"):
                parse_folder(folder)
        
        st.divider()
        
        # Data status
        df = get_df()
        if df is not None:
            st.success(f"✅ {len(df)} molecules loaded")
            
            # Quick filters
            st.subheader("🔍 Quick Filters")
            
            # Molecule filter
            mol_ids = df["molecule_id"].dropna().unique().tolist()
            selected_mols = st.multiselect(
                "Filter Molecules",
                mol_ids,
                default=mol_ids[:min(5, len(mol_ids))],
                key="global_mol_filter"
            )
            st.session_state.selected_molecules = selected_mols
            
            # State filter
            if "optimized_state" in df.columns:
                states = df["optimized_state"].dropna().unique().tolist()
                selected_states = st.multiselect(
                    "Filter States",
                    states,
                    default=states,
                    key="global_state_filter"
                )
                st.session_state.selected_states = selected_states
    
    # Main content
    df = get_df()
    if df is not None:
        show_main_content(df)
    else:
        show_welcome()


def show_welcome():
    """Welcome screen."""
    st.title("🔬 ORCA Quantum Chemistry Visualization")
    st.markdown("**Interactive analysis and visualization of ORCA output files**")
    
    st.info("👈 Upload files or download sample data to get started.")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown("""
        ### 📊 Visualizations
        - Energy diagrams
        - IR/Raman/UV-Vis spectra
        - Orbital energy levels
        - 2D/3D molecular viewer
        """)
    
    with col2:
        st.markdown("""
        ### 🔧 Features
        - Interactive customization
        - Multi-molecule comparison
        - Adaptive data tables
        - Export to CSV/Excel/HTML
        """)
    
    with col3:
        st.markdown("""
        ### 🔗 Reaction Editor
        - Node-based visualization
        - Drag-and-drop molecules
        - Connect reaction pathways
        - Click nodes for details
        """)


def show_main_content(df):
    """Main content with tabs."""
    
    # Apply global filters
    filtered_df = df.copy()
    if hasattr(st.session_state, 'selected_molecules') and st.session_state.selected_molecules:
        filtered_df = filtered_df[filtered_df["molecule_id"].isin(st.session_state.selected_molecules)]
    
    # Tabs
    tabs = st.tabs([
        "📊 Dashboard",
        "🧬 Molecules", 
        "📈 Spectra",
        "⚡ Energy",
        "🔮 Orbitals",
        "🔗 Reactions",
        "📤 Export"
    ])
    
    with tabs[0]:
        show_dashboard(filtered_df)
    
    with tabs[1]:
        show_molecules_tab(filtered_df)
    
    with tabs[2]:
        show_spectra_tab(filtered_df)
    
    with tabs[3]:
        show_energy_tab(filtered_df)
    
    with tabs[4]:
        show_orbitals_tab(filtered_df)
    
    with tabs[5]:
        show_reactions_tab(filtered_df)
    
    with tabs[6]:
        show_export_tab(filtered_df)


def show_dashboard(df):
    """Dashboard overview."""
    st.header("📊 Dashboard")
    
    # Metrics
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Molecules", len(df))
    
    with col2:
        has_energy = df["gibbs_Eh"].notna().sum()
        st.metric("With Energy", has_energy)
    
    with col3:
        has_ir = df["ir"].apply(lambda x: x is not None and hasattr(x, '__len__') and len(x) > 0).sum()
        st.metric("With IR", has_ir)
    
    with col4:
        has_tddft = df["has_tddft"].sum() if "has_tddft" in df.columns else 0
        st.metric("With TDDFT", has_tddft)
    
    st.divider()
    
    # Quick data table
    st.subheader("📋 Data Overview")
    scalar_cols = ["molecule_id", "method_id", "functional", "basis_set", 
                   "gibbs_Eh", "single_point_Eh", "homo_energy", "lumo_energy",
                   "optimized_state", "calc_class"]
    available = [c for c in scalar_cols if c in df.columns]
    st.dataframe(df[available], use_container_width=True, height=400)


def show_molecules_tab(df):
    """Molecule viewer tab."""
    from components.molecule_viewer import render_molecule_viewer
    render_molecule_viewer(df)


def show_spectra_tab(df):
    """Spectra visualization tab."""
    from components.spectra_viz import render_spectra_viz
    render_spectra_viz(df)


def show_energy_tab(df):
    """Energy comparison tab."""
    from components.energy_viz import render_energy_viz
    render_energy_viz(df)


def show_orbitals_tab(df):
    """Orbital visualization tab."""
    from components.orbital_viz import render_orbital_viz
    render_orbital_viz(df)


def show_reactions_tab(df):
    """Reaction node editor tab."""
    from components.node_editor import render_node_editor
    render_node_editor(df)


def show_export_tab(df):
    """Export tab."""
    from components.export_panel import render_export_panel
    render_export_panel(df)


# === Data Loading Functions ===

def parse_uploaded_files(files):
    """Parse uploaded files."""
    from src.parser.factory import ParserFactory
    
    factory = ParserFactory()
    results = []
    
    progress = st.progress(0, "Parsing files...")
    for i, file in enumerate(files):
        try:
            content = file.read().decode('utf-8', errors='ignore')
            result = factory.parse_text(content, file.name)
            results.append(result.to_dict())
        except Exception as e:
            st.warning(f"Failed: {file.name} - {e}")
        progress.progress((i + 1) / len(files))
    
    if results:
        import pandas as pd
        set_df(pd.DataFrame(results))
        st.success(f"✅ Parsed {len(results)} files")
        st.rerun()


def download_hf_data():
    """Download from HuggingFace and auto-parse."""
    try:
        from huggingface_hub import snapshot_download
        
        st.info("📥 Step 1/2: Downloading from HuggingFace...")
        with st.spinner("Downloading..."):
            snapshot_download(
                repo_id="JauharMz/Orca",
                repo_type="dataset",
                local_dir="./test_data_hf",
                local_dir_use_symlinks=False
            )
        st.success("✅ Downloaded to ./test_data_hf")
        
        # Auto-parse immediately
        st.info("🔄 Step 2/2: Parsing files...")
        parse_folder("./test_data_hf")
        
    except ImportError:
        st.error("❌ huggingface_hub not installed. Run: pip install huggingface_hub")
    except Exception as e:
        st.error(f"❌ Download failed: {e}")


def parse_folder(folder):
    """Parse folder of .out files."""
    from src.parser.batch import BatchParser
    import pandas as pd
    
    pattern = f"{folder}/**/*.out"
    
    with st.spinner(f"Parsing {pattern}..."):
        batch = BatchParser(pattern)
        df = batch.parse_all(verbose=False)
    
    if df is not None and len(df) > 0:
        set_df(df)
        st.success(f"✅ Parsed {len(df)} files")
        st.rerun()
    else:
        st.warning("No files found or parsing failed")


if __name__ == "__main__":
    main()
