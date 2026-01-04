"""
ORCA Visualization Platform v2 - Main Entry

Run: streamlit run streamlit_app/app.py
"""

import streamlit as st
import sys
import os
import glob
import re
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
    .element-container { margin-bottom: 0.3rem; }
</style>
""", unsafe_allow_html=True)

# Import utilities
from utils.session_state import init_session_state, get_df, set_df

# Initialize
init_session_state()


def main():
    """Main app entry."""
    
    # Sidebar - Data Loading only
    with st.sidebar:
        st.title("🔬 ORCA Platform")
        st.divider()
        
        st.header("📁 Data Loading")
        
        data_source = st.radio("Source", ["HuggingFace", "Local Folder", "Upload"], horizontal=True)
        
        if data_source == "HuggingFace":
            data_path = Path("./test_data_hf")
            if data_path.exists() and list(data_path.rglob("*.out")):
                file_count = len(list(data_path.rglob("*.out")))
                st.info(f"📁 {file_count} files in test_data_hf")
                if st.button("🚀 Parse Data", type="primary"):
                    parse_folder_per_file("./test_data_hf")
            else:
                if st.button("📥 Download & Parse", type="primary"):
                    download_hf_data()
        
        elif data_source == "Local Folder":
            folder = st.text_input("Folder", "./test_data_hf")
            if st.button("📂 Parse", type="primary"):
                parse_folder_per_file(folder)
        
        elif data_source == "Upload":
            files = st.file_uploader("Upload .out files", type=["out"], accept_multiple_files=True)
            if files and st.button("🚀 Parse", type="primary"):
                parse_uploaded_files(files)
        
        st.divider()
        
        # Status only
        df = get_df()
        if df is not None:
            st.success(f"✅ {len(df)} records loaded")
            
            # Show state distribution
            if "optimized_state" in df.columns:
                states = df["optimized_state"].value_counts()
                st.caption("States: " + ", ".join([f"{s}={c}" for s, c in states.items()]))
    
    # Main content
    df = get_df()
    if df is not None:
        show_main_content(df)
    else:
        show_welcome()


def show_welcome():
    """Welcome screen."""
    st.title("🔬 ORCA Quantum Chemistry Visualization")
    st.markdown("**Interactive analysis of ORCA output files**")
    st.info("👈 Load data from the sidebar to get started.")
    
    col1, col2, col3 = st.columns(3)
    with col1:
        st.markdown("### 📊 Visualizations\n- Energy diagrams\n- IR/Raman/UV-Vis\n- Orbital levels")
    with col2:
        st.markdown("### 🔧 Features\n- Multi-molecule compare\n- Data export\n- Customization")
    with col3:
        st.markdown("### 🔗 Reactions\n- Node-based editor\n- Energy pathways")


def show_main_content(df):
    """Main content - no top-level filter, filters are per-visualization."""
    
    st.header("🔬 ORCA Visualization")
    
    # Quick stats
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Total Records", len(df))
    with col2:
        unique_mols = df["molecule_id"].nunique() if "molecule_id" in df.columns else 0
        st.metric("Molecules", unique_mols)
    with col3:
        if "optimized_state" in df.columns:
            unique_states = df["optimized_state"].nunique()
            st.metric("States", unique_states)
    with col4:
        has_coords = df["cart_coords"].apply(lambda x: x is not None and hasattr(x, '__len__') and len(x) > 0).sum() if "cart_coords" in df.columns else 0
        st.metric("With Coords", has_coords)
    
    # Tabs - each tab has its own filter
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
    """Data overview with filtering."""
    st.subheader("📊 Data Overview")
    
    # Filters
    col1, col2 = st.columns(2)
    with col1:
        mol_ids = sorted(df["molecule_id"].dropna().unique().tolist()) if "molecule_id" in df.columns else []
        selected_mols = st.multiselect("Filter by Molecule", mol_ids, default=mol_ids[:min(10, len(mol_ids))], key="data_mol_filter")
    with col2:
        states = sorted(df["optimized_state"].dropna().unique().tolist()) if "optimized_state" in df.columns else []
        selected_states = st.multiselect("Filter by State", states, default=states, key="data_state_filter")
    
    # Apply filter
    filtered = df.copy()
    if selected_mols:
        filtered = filtered[filtered["molecule_id"].isin(selected_mols)]
    if selected_states and "optimized_state" in filtered.columns:
        filtered = filtered[filtered["optimized_state"].isin(selected_states)]
    
    st.caption(f"Showing {len(filtered)} of {len(df)} records")
    
    # Show table
    display_cols = ["molecule_id", "optimized_state", "method_id", "functional", "basis_set", 
                    "gibbs_Eh", "single_point_Eh", "homo_energy", "lumo_energy"]
    available = [c for c in display_cols if c in filtered.columns]
    st.dataframe(filtered[available], use_container_width=True, height=500)


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


# === Data Loading Functions ===

def parse_uploaded_files(files):
    """Parse uploaded files."""
    from src.parser.factory import ParserFactory
    import pandas as pd
    
    factory = ParserFactory()
    results = []
    
    progress = st.progress(0, "Parsing...")
    for i, file in enumerate(files):
        try:
            content = file.read().decode('utf-8', errors='ignore')
            result = factory.parse_text(content, file.name)
            row = result.to_dict()
            row["molecule_id"] = extract_mol_id(file.name)
            row["optimized_state"] = extract_state(file.name)
            results.append(row)
        except Exception:
            pass
        progress.progress((i + 1) / len(files))
    
    if results:
        set_df(pd.DataFrame(results))
        st.success(f"✅ Parsed {len(results)} files")
        st.rerun()


def download_hf_data():
    """Download from HuggingFace."""
    try:
        from huggingface_hub import snapshot_download
        
        st.info("📥 Downloading...")
        with st.spinner("Downloading..."):
            snapshot_download(
                repo_id="JauharMz/Orca",
                repo_type="dataset",
                local_dir="./test_data_hf",
                local_dir_use_symlinks=False
            )
        st.success("✅ Downloaded")
        parse_folder_per_file("./test_data_hf")
        
    except ImportError:
        st.error("Install: pip install huggingface_hub")
    except Exception as e:
        st.error(f"Failed: {e}")


def parse_folder_per_file(folder: str):
    """Parse each file as a separate record with state from folder name."""
    from src.parser.factory import ParserFactory
    import pandas as pd
    
    pattern = f"{folder}/**/*.out"
    files = sorted(glob.glob(pattern, recursive=True))
    
    if not files:
        st.warning(f"No .out files in {folder}")
        return
    
    st.info(f"🔄 Parsing {len(files)} files...")
    
    factory = ParserFactory()
    results = []
    
    progress = st.progress(0, "Parsing...")
    for i, filepath in enumerate(files):
        try:
            result = factory.parse(filepath)
            row = result.to_dict()
            
            # Extract molecule_id and state from path
            mol_id, state = extract_molecule_and_state(filepath)
            row["molecule_id"] = mol_id
            row["optimized_state"] = state
            row["_filepath"] = filepath
            
            results.append(row)
        except Exception:
            pass
        progress.progress((i + 1) / len(files))
    
    if not results:
        st.warning("No files parsed")
        return
    
    df = pd.DataFrame(results)
    set_df(df)
    
    # Stats
    if "optimized_state" in df.columns:
        states = df["optimized_state"].value_counts()
        state_str = ", ".join([f"{s}={c}" for s, c in states.items()])
        st.success(f"✅ Parsed {len(df)} records: {state_str}")
    else:
        st.success(f"✅ Parsed {len(df)} records")
    st.rerun()


def extract_molecule_and_state(filepath: str) -> tuple:
    """
    Extract molecule ID and state from filepath.
    
    Examples:
        p1as0/file.out -> (p1a, S0)
        p1xs0p/file.out -> (p1x, S0-SP)
        p1xvg/file.out -> (p1x, VG)
    """
    parts = Path(filepath).parts
    
    # Pattern to match state suffix
    state_pattern = re.compile(r'(s0p|s0|s1|t1|vg|ah|ahas)$', re.I)
    
    # Try parent folder first
    for part in reversed(parts[:-1]):  # All except filename
        match = state_pattern.search(part.lower())
        if match:
            state_code = match.group(1).upper()
            mol_id = state_pattern.sub('', part).strip('_-')
            
            # Map to readable names
            state_map = {
                "S0P": "S0-SP",
                "S0": "S0", 
                "S1": "S1",
                "T1": "T1",
                "VG": "VG",
                "AH": "AH",
                "AHAS": "AHAS"
            }
            state = state_map.get(state_code, state_code)
            return mol_id, state
    
    # Fallback: use filename
    filename = os.path.splitext(os.path.basename(filepath))[0]
    match = state_pattern.search(filename.lower())
    if match:
        state_code = match.group(1).upper()
        mol_id = state_pattern.sub('', filename).strip('_-')
        return mol_id, state_code
    
    return filename, "unknown"


def extract_mol_id(filename: str) -> str:
    """Extract molecule ID from filename."""
    base = os.path.splitext(os.path.basename(filename))[0]
    mol_id = re.sub(r'(_?s0p?|_?s0|_?s1|_?t1|_?vg|_?opt)$', '', base, flags=re.I)
    return mol_id.strip('_-')


def extract_state(filename: str) -> str:
    """Extract state from filename."""
    match = re.search(r'(s0p|s0|s1|t1|vg)$', filename.lower())
    if match:
        return match.group(1).upper()
    return "unknown"


if __name__ == "__main__":
    main()
