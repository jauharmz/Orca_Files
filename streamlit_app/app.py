"""
ORCA Visualization Platform - Main Entry

Uses ParserOrchestrator for data loading.
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

# Custom CSS
st.markdown("""
<style>
    .block-container { padding-top: 1rem; padding-bottom: 1rem; }
    .stTabs [data-baseweb="tab-list"] { gap: 8px; }
</style>
""", unsafe_allow_html=True)

# Import orchestrator and utilities
from src.orchestrator.parser_orchestrator import ParserOrchestrator
from utils.session_state import init_session_state, get_df, set_df

# Initialize
init_session_state()


def get_orchestrator() -> ParserOrchestrator:
    """Get or create orchestrator instance."""
    if "orchestrator" not in st.session_state:
        st.session_state.orchestrator = ParserOrchestrator()
    return st.session_state.orchestrator


def main():
    """Main app entry."""
    
    # Sidebar - Data Loading
    with st.sidebar:
        st.title("🔬 ORCA Platform")
        st.divider()
        
        st.header("📁 Data Loading")
        
        data_source = st.radio("Source", ["Sample Data", "Local Folder"], horizontal=True)
        
        orchestrator = get_orchestrator()
        
        if data_source == "Sample Data":
            data_path = Path("./test_data_hf")
            if data_path.exists() and list(data_path.rglob("*.out")):
                file_count = len(list(data_path.rglob("*.out")))
                st.info(f"📁 {file_count} files found")
                if st.button("🚀 Parse Data", type="primary"):
                    with st.spinner("Parsing..."):
                        df = orchestrator.parse_folder("./test_data_hf")
                        set_df(df)
                    st.rerun()
            else:
                if st.button("📥 Download Sample", type="primary"):
                    download_sample_data()
        
        elif data_source == "Local Folder":
            folder = st.text_input("Folder Path", "./test_data_hf")
            if st.button("📂 Parse Folder", type="primary"):
                with st.spinner("Parsing..."):
                    df = orchestrator.parse_folder(folder)
                    set_df(df)
                st.rerun()
        
        st.divider()
        
        # Show data summary
        df = get_df()
        if df is not None and len(df) > 0:
            summary = orchestrator.get_summary()
            st.success(f"✅ {summary['total']} records")
            st.caption(f"Molecules: {summary.get('molecules', 0)}")
            if "states" in summary:
                state_str = ", ".join([f"{s}={c}" for s, c in summary["states"].items()])
                st.caption(f"States: {state_str}")
    
    # Main content
    df = get_df()
    if df is not None and len(df) > 0:
        show_main_content(df)
    else:
        show_welcome()


def show_welcome():
    """Welcome screen."""
    st.title("🔬 ORCA Quantum Chemistry Visualization")
    st.markdown("**Interactive analysis of ORCA output files**")
    st.info("👈 Load data from the sidebar to get started.")


def show_main_content(df):
    """Main content with per-tab molecule/state filtering."""
    
    st.header("🔬 ORCA Visualization")
    
    # Quick stats
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        st.metric("Records", len(df))
    with c2:
        st.metric("Molecules", df["molecule_id"].nunique() if "molecule_id" in df.columns else 0)
    with c3:
        if "optimized_state" in df.columns:
            st.metric("States", df["optimized_state"].nunique())
    with c4:
        has_coords = df["cart_coords"].apply(lambda x: x is not None and hasattr(x, '__len__') and len(x) > 0).sum() if "cart_coords" in df.columns else 0
        st.metric("With Coords", has_coords)
    
    # Tabs
    tabs = st.tabs([
        "📊 Data",
        "🧬 Molecules",
        "📈 Spectra",
        "⚡ Energy",
        "🔮 Orbitals",
        "🔗 Reactions",
        "📤 Export"
    ])
    
    with tabs[0]:
        show_data_tab(df)
    with tabs[1]:
        show_molecules_tab(df)
    with tabs[2]:
        show_spectra_tab(df)
    with tabs[3]:
        show_energy_tab(df)
    with tabs[4]:
        show_orbitals_tab(df)
    with tabs[5]:
        show_reactions_tab(df)
    with tabs[6]:
        show_export_tab(df)


def show_data_tab(df):
    """Data overview with molecule+state filtering."""
    st.subheader("📊 Data Overview")
    
    # Dual filters
    c1, c2 = st.columns(2)
    with c1:
        mol_ids = sorted(df["molecule_id"].dropna().unique().tolist()) if "molecule_id" in df.columns else []
        sel_mols = st.multiselect("Molecules", mol_ids, default=mol_ids[:min(10, len(mol_ids))], key="data_mols")
    with c2:
        states = sorted(df["optimized_state"].dropna().unique().tolist()) if "optimized_state" in df.columns else []
        sel_states = st.multiselect("States", states, default=states, key="data_states")
    
    # Filter
    filtered = df.copy()
    if sel_mols:
        filtered = filtered[filtered["molecule_id"].isin(sel_mols)]
    if sel_states and "optimized_state" in filtered.columns:
        filtered = filtered[filtered["optimized_state"].isin(sel_states)]
    
    st.caption(f"Showing {len(filtered)} of {len(df)} records")
    
    # Display columns
    cols = ["molecule_id", "optimized_state", "functional", "basis_set", 
            "gibbs_Eh", "single_point_Eh", "homo_energy", "lumo_energy", "calc_class"]
    disp = [c for c in cols if c in filtered.columns]
    st.dataframe(filtered[disp], use_container_width=True, height=500)


def show_molecules_tab(df):
    """Molecule viewer."""
    from components.molecule_viewer import render_molecule_viewer
    render_molecule_viewer(df)


def show_spectra_tab(df):
    """Spectra visualization."""
    from components.spectra_viz import render_spectra_viz
    render_spectra_viz(df)


def show_energy_tab(df):
    """Energy comparison."""
    from components.energy_viz import render_energy_viz
    render_energy_viz(df)


def show_orbitals_tab(df):
    """Orbital visualization."""
    from components.orbital_viz import render_orbital_viz
    render_orbital_viz(df)


def show_reactions_tab(df):
    """Reaction node editor."""
    from components.node_editor import render_node_editor
    render_node_editor(df)


def show_export_tab(df):
    """Export panel."""
    from components.export_panel import render_export_panel
    render_export_panel(df)


def download_sample_data():
    """Download sample data from HuggingFace."""
    try:
        from huggingface_hub import snapshot_download
        
        st.info("📥 Downloading from HuggingFace...")
        with st.spinner("Downloading..."):
            snapshot_download(
                repo_id="JauharMz/Orca",
                repo_type="dataset",
                local_dir="./test_data_hf",
                local_dir_use_symlinks=False
            )
        st.success("✅ Downloaded to ./test_data_hf")
        st.rerun()
        
    except ImportError:
        st.error("Install: pip install huggingface_hub")
    except Exception as e:
        st.error(f"Failed: {e}")


if __name__ == "__main__":
    main()
