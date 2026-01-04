"""
ORCA Visualization Platform v2

No top-level filters - each tab has its own filter.
Run: streamlit run streamlit_app/app.py
"""

import streamlit as st
import sys
import os
import glob
import re
from pathlib import Path
import pandas as pd

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

# Import utilities
from utils.session_state import init_session_state, get_df, set_df

# Initialize
init_session_state()


def main():
    """Main app entry."""
    
    # Sidebar - Data Loading
    with st.sidebar:
        st.title("🔬 ORCA Platform")
        st.divider()
        
        st.header("📁 Data Loading")
        
        data_source = st.radio(
            "Source",
            ["HuggingFace Sample", "Local Folder", "Upload Files"],
            horizontal=True
        )
        
        if data_source == "HuggingFace Sample":
            data_path = Path("./test_data_hf")
            if data_path.exists():
                file_count = len(list(data_path.rglob("*.out")))
                st.info(f"📁 Sample data ({file_count} files)")
                if st.button("🚀 Parse Data", type="primary"):
                    parse_folder("./test_data_hf")
            else:
                if st.button("📥 Download & Parse", type="primary"):
                    download_hf_data()
        
        elif data_source == "Local Folder":
            folder_path = st.text_input("Folder Path", value="./test_data_hf")
            if st.button("📂 Parse Folder", type="primary"):
                parse_folder(folder_path)
        
        elif data_source == "Upload Files":
            uploaded_files = st.file_uploader(
                "Upload .out files",
                type=["out"],
                accept_multiple_files=True
            )
            if uploaded_files and st.button("🚀 Parse Uploaded", type="primary"):
                parse_uploaded_files(uploaded_files)
        
        st.divider()
        
        # Data Status
        df = get_df()
        if df is not None and len(df) > 0:
            st.success(f"✅ {len(df)} records loaded")
            
            if "optimized_state" in df.columns:
                states = df["optimized_state"].value_counts()
                st.caption("States: " + ", ".join([f"{s}={c}" for s, c in states.items()]))
            
            if "molecule_id" in df.columns:
                st.caption(f"Molecules: {df['molecule_id'].nunique()}")
    
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
    """Main content - NO top-level filters, each tab handles its own."""
    
    st.header("🔬 ORCA Visualization Platform")
    
    # Quick stats only (no filters)
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        st.metric("Records", len(df))
    with c2:
        st.metric("Molecules", df["molecule_id"].nunique() if "molecule_id" in df.columns else 0)
    with c3:
        if "optimized_state" in df.columns:
            st.metric("States", df["optimized_state"].nunique())
    with c4:
        if "gibbs_Eh" in df.columns:
            st.metric("With Energy", df["gibbs_Eh"].notna().sum())
    
    # Main tabs - each handles its own filtering
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
    """Data table."""
    from components.data_table import render_data_table
    render_data_table(df)


def show_molecules_tab(df):
    """Molecule viewer."""
    from components.molecule_viewer import render_molecule_viewer
    render_molecule_viewer(df)


def show_spectra_tab(df):
    """Spectra visualization."""
    from components.spectra_viz import render_spectra_viz
    render_spectra_viz(df)


def show_energy_tab(df):
    """Energy visualization."""
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
    
    factory = ParserFactory()
    results = []
    
    progress = st.progress(0, "Parsing files...")
    for i, file in enumerate(files):
        try:
            content = file.read().decode('utf-8', errors='ignore')
            result = factory.parse_text(content, file.name)
            row = result.to_dict()
            
            mol_id, state = extract_molecule_and_state(file.name)
            row["molecule_id"] = mol_id
            row["optimized_state"] = state
            
            results.append(row)
        except Exception as e:
            st.warning(f"Failed: {file.name} - {e}")
        progress.progress((i + 1) / len(files))
    
    if results:
        set_df(pd.DataFrame(results))
        st.success(f"✅ Parsed {len(results)} files")
        st.rerun()


def download_hf_data():
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
        st.success("✅ Downloaded")
        parse_folder("./test_data_hf")
        
    except ImportError:
        st.error("Install: pip install huggingface_hub")
    except Exception as e:
        st.error(f"Failed: {e}")


def parse_folder(folder: str):
    """Parse folder using BatchParser.parse_all_detailed()."""
    from src.parser.batch import BatchParser
    
    pattern = f"{folder}/**/*.out"
    files = sorted(glob.glob(pattern, recursive=True))
    
    if not files:
        st.warning(f"No .out files found in {folder}")
        return
    
    st.info(f"🔄 Parsing {len(files)} files...")
    
    # Show progress
    progress = st.progress(0, "Parsing...")
    
    batch = BatchParser(pattern)
    
    # Use detailed mode for state extraction
    if hasattr(batch, 'parse_all_detailed'):
        df = batch.parse_all_detailed(verbose=False)
    else:
        df = batch.parse_all(verbose=False)
    
    progress.progress(100)
    
    if df is not None and len(df) > 0:
        set_df(df)
        
        # Show summary
        if "optimized_state" in df.columns:
            states = df["optimized_state"].value_counts()
            st.success(f"✅ Parsed {len(df)} records: " + ", ".join([f"{s}={c}" for s, c in states.items()]))
        else:
            st.success(f"✅ Parsed {len(df)} records")
        st.rerun()
    else:
        st.warning("No files parsed successfully")


def extract_molecule_and_state(filepath: str) -> tuple:
    """Extract molecule ID and state from filepath."""
    parts = Path(filepath).parts
    state_pattern = re.compile(r'(s0p|s0|s1|t1|vg|ah|ahas)(?:[^a-z]|$)', re.I)
    
    for part in reversed(parts[:-1]):
        match = state_pattern.search(part.lower())
        if match:
            state_code = match.group(1).upper()
            mol_id = state_pattern.sub('', part).strip('_-')
            state_map = {"S0P": "S0-SP", "S0": "S0", "S1": "S1", "T1": "T1", "VG": "VG", "AH": "AH", "AHAS": "AHAS"}
            return mol_id, state_map.get(state_code, state_code)
    
    filename = os.path.splitext(os.path.basename(filepath))[0]
    match = state_pattern.search(filename.lower())
    if match:
        state_code = match.group(1).upper()
        mol_id = state_pattern.sub('', filename).strip('_-')
        return mol_id, state_code
    
    return filename, "unknown"


if __name__ == "__main__":
    main()
