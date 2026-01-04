"""
ORCA Visualization Platform v2 - Complete Implementation

Follows architecture from README.md:
- Parser Layer → Analysis Layer → Data Layer → Visualization Layer → Export Layer

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

# Custom CSS for compact layout
st.markdown("""
<style>
    .block-container { padding-top: 1rem; padding-bottom: 1rem; }
    .element-container { margin-bottom: 0.3rem; }
    .stTabs [data-baseweb="tab-list"] { gap: 8px; }
    div[data-testid="stMetric"] { background: #f8f9fa; padding: 10px; border-radius: 8px; }
</style>
""", unsafe_allow_html=True)

# Import utilities
from utils.session_state import init_session_state, get_df, set_df

# Initialize
init_session_state()


def main():
    """Main app entry."""
    
    # Sidebar
    with st.sidebar:
        st.title("🔬 ORCA Platform")
        st.divider()
        
        # Data Source Selection
        st.header("📁 Data Loading")
        
        data_source = st.radio(
            "Select Source",
            ["HuggingFace Sample", "Local Folder", "Upload Files"],
            horizontal=True
        )
        
        if data_source == "HuggingFace Sample":
            data_path = Path("./test_data_hf")
            if data_path.exists():
                file_count = len(list(data_path.rglob("*.out")))
                st.info(f"📁 Sample data exists ({file_count} files)")
                col1, col2 = st.columns(2)
                with col1:
                    if st.button("🚀 Parse Existing", type="primary"):
                        parse_folder("./test_data_hf")
                with col2:
                    if st.button("🔄 Re-download"):
                        download_hf_data()
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
            if uploaded_files:
                if st.button("🚀 Parse Uploaded", type="primary"):
                    parse_uploaded_files(uploaded_files)
        
        st.divider()
        
        # Data Status
        df = get_df()
        if df is not None and len(df) > 0:
            st.success(f"✅ {len(df)} molecules loaded")
            
            # State distribution
            if "optimized_state" in df.columns:
                states = df["optimized_state"].value_counts()
                st.caption("📊 States: " + ", ".join([f"{s}={c}" for s, c in states.items()]))
            
            # Quick stats
            with st.expander("📈 Quick Stats", expanded=False):
                st.write(f"Unique molecules: {df['molecule_id'].nunique()}")
                for col in ["gibbs_Eh", "homo_energy", "cart_coords"]:
                    if col in df.columns:
                        count = df[col].notna().sum()
                        st.write(f"{col}: {count}/{len(df)}")
        
        st.divider()
        
        # Settings
        with st.expander("⚙️ Settings", expanded=False):
            st.checkbox("Debug mode", key="debug_mode")
            st.checkbox("Show raw data", key="show_raw")
    
    # Main content
    df = get_df()
    if df is not None and len(df) > 0:
        show_main_content(df)
    else:
        show_welcome()


def show_welcome():
    """Welcome screen when no data loaded."""
    st.title("🔬 ORCA Quantum Chemistry Visualization")
    
    st.markdown("""
    **An interactive platform for analyzing ORCA output files.**
    
    ### Features:
    - 📊 **Data Parsing**: Modular parsers for energy, geometry, orbitals, spectra
    - 🧬 **3D Visualization**: Interactive molecular structures with py3Dmol
    - 📈 **Spectroscopy**: IR, Raman, UV-Vis with smoothing and comparison
    - ⚡ **Energy Analysis**: HOMO-LUMO diagrams, relative energies
    - 🔗 **Pathway Detection**: Reaction pathways with node editor
    - 📤 **Export**: CSV, Excel, JSON, HTML reports
    """)
    
    st.info("👈 **Load data from the sidebar to get started.**")
    
    # Quick actions
    col1, col2, col3 = st.columns(3)
    with col1:
        st.markdown("### 📥 Sample Data")
        st.markdown("Download from HuggingFace")
    with col2:
        st.markdown("### 📂 Local Files")
        st.markdown("Parse .out files from folder")
    with col3:
        st.markdown("### 📤 Upload")
        st.markdown("Upload files directly")


def show_main_content(df):
    """Main content with all visualizations."""
    
    st.header("🔬 ORCA Visualization Platform")
    
    # Global filters at top
    with st.container():
        col1, col2, col3, col4 = st.columns([3, 3, 1, 1])
        
        # Molecule filter
        with col1:
            mol_ids = sorted(df["molecule_id"].dropna().unique().tolist())
            
            # Initialize or validate session state
            if "selected_molecules" not in st.session_state:
                st.session_state.selected_molecules = mol_ids[:min(10, len(mol_ids))]
            else:
                st.session_state.selected_molecules = [
                    m for m in st.session_state.selected_molecules if m in mol_ids
                ]
            
            valid_defaults = [m for m in st.session_state.selected_molecules if m in mol_ids]
            selected_mols = st.multiselect(
                "🔍 Filter Molecules",
                mol_ids,
                default=valid_defaults if valid_defaults else mol_ids[:min(5, len(mol_ids))],
                key="mol_filter_widget"
            )
            st.session_state.selected_molecules = selected_mols
        
        # State filter
        with col2:
            if "optimized_state" in df.columns:
                states = sorted(df["optimized_state"].dropna().unique().tolist())
                if states:
                    if "selected_states" not in st.session_state:
                        st.session_state.selected_states = states
                    valid_states = [s for s in st.session_state.selected_states if s in states]
                    selected_states = st.multiselect(
                        "📊 Filter States", 
                        states, 
                        default=valid_states if valid_states else states,
                        key="state_filter_widget"
                    )
                    st.session_state.selected_states = selected_states
                else:
                    selected_states = []
            else:
                selected_states = []
        
        # Quick filter buttons
        with col3:
            if st.button("✅ All", key="select_all"):
                st.session_state.selected_molecules = mol_ids
                st.rerun()
        with col4:
            if st.button("❌ Clear", key="clear_all"):
                st.session_state.selected_molecules = []
                st.rerun()
    
    # Apply filters
    filtered_df = df.copy()
    if selected_mols:
        filtered_df = filtered_df[filtered_df["molecule_id"].isin(selected_mols)]
    if selected_states and "optimized_state" in filtered_df.columns:
        filtered_df = filtered_df[filtered_df["optimized_state"].isin(selected_states)]
    
    st.caption(f"Showing {len(filtered_df)} of {len(df)} records")
    
    # Stats row
    c1, c2, c3, c4, c5 = st.columns(5)
    with c1:
        st.metric("Records", len(filtered_df))
    with c2:
        st.metric("Molecules", filtered_df["molecule_id"].nunique())
    with c3:
        if "optimized_state" in filtered_df.columns:
            st.metric("States", filtered_df["optimized_state"].nunique())
    with c4:
        if "gibbs_Eh" in filtered_df.columns:
            st.metric("With Energy", filtered_df["gibbs_Eh"].notna().sum())
    with c5:
        if "cart_coords" in filtered_df.columns:
            has_coords = filtered_df["cart_coords"].apply(
                lambda x: x is not None and hasattr(x, '__len__') and len(x) > 0
            ).sum()
            st.metric("With Coords", has_coords)
    
    # Main tabs
    tabs = st.tabs([
        "📊 Data Table",
        "🧬 Molecules",
        "📈 Spectra",
        "⚡ Energy",
        "🔮 Orbitals",
        "🔗 Reactions",
        "📤 Export"
    ])
    
    with tabs[0]:
        show_data_table_tab(filtered_df)
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


def show_data_table_tab(df):
    """Data table with all columns."""
    from components.data_table import render_data_table
    render_data_table(df)


def show_molecules_tab(df):
    """Molecule 2D/3D viewer."""
    from components.molecule_viewer import render_molecule_viewer
    render_molecule_viewer(df)


def show_spectra_tab(df):
    """IR, Raman, UV-Vis spectra."""
    from components.spectra_viz import render_spectra_viz
    render_spectra_viz(df)


def show_energy_tab(df):
    """Energy comparison and diagrams."""
    from components.energy_viz import render_energy_viz
    render_energy_viz(df)


def show_orbitals_tab(df):
    """Orbital energy levels."""
    from components.orbital_viz import render_orbital_viz
    render_orbital_viz(df)


def show_reactions_tab(df):
    """Reaction pathway node editor."""
    from components.node_editor import render_node_editor
    render_node_editor(df)


def show_export_tab(df):
    """Export panel."""
    from components.export_panel import render_export_panel
    render_export_panel(df)


# === Data Loading Functions ===

def parse_uploaded_files(files):
    """Parse uploaded files - one record per file with state extraction."""
    from src.parser.factory import ParserFactory
    
    factory = ParserFactory()
    results = []
    
    progress = st.progress(0, "Parsing files...")
    for i, file in enumerate(files):
        try:
            content = file.read().decode('utf-8', errors='ignore')
            result = factory.parse_text(content, file.name)
            row = result.to_dict()
            
            # Extract molecule_id and state from filename
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
    """Download sample data from HuggingFace and parse."""
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
        
        st.info("🔄 Step 2/2: Parsing files...")
        parse_folder("./test_data_hf")
        
    except ImportError:
        st.error("❌ huggingface_hub not installed. Run: pip install huggingface_hub")
    except Exception as e:
        st.error(f"❌ Download failed: {e}")


def parse_folder(folder: str):
    """Parse folder - uses detailed mode for state extraction."""
    from src.parser.batch import BatchParser
    
    pattern = f"{folder}/**/*.out"
    files = sorted(glob.glob(pattern, recursive=True))
    
    if not files:
        st.warning(f"No .out files found in {folder}")
        return
    
    st.info(f"🔄 Parsing {len(files)} files...")
    
    batch = BatchParser(pattern)
    
    # Use detailed mode if available, otherwise fall back to standard
    if hasattr(batch, 'parse_all_detailed'):
        df = batch.parse_all_detailed(verbose=True)
    else:
        df = batch.parse_all(verbose=True)
    
    if df is not None and len(df) > 0:
        set_df(df)
        
        # Show state distribution
        if "optimized_state" in df.columns:
            states = df["optimized_state"].value_counts()
            state_str = ", ".join([f"{s}={c}" for s, c in states.items()])
            st.success(f"✅ Parsed {len(df)} records: {state_str}")
        else:
            st.success(f"✅ Parsed {len(df)} molecules")
        st.rerun()
    else:
        st.warning("No files parsed successfully")


def extract_molecule_and_state(filepath: str) -> tuple:
    """Extract molecule ID and state from filepath."""
    parts = Path(filepath).parts
    
    # Pattern to match state suffix in folder or filename
    state_pattern = re.compile(r'(s0p|s0|s1|t1|vg|ah|ahas)(?:[^a-z]|$)', re.I)
    
    # Try folder names first
    for part in reversed(parts[:-1]):
        match = state_pattern.search(part.lower())
        if match:
            state_code = match.group(1).upper()
            mol_id = state_pattern.sub('', part).strip('_-')
            
            state_map = {
                "S0P": "S0-SP", "S0": "S0", "S1": "S1", "T1": "T1",
                "VG": "VG", "AH": "AH", "AHAS": "AHAS"
            }
            return mol_id, state_map.get(state_code, state_code)
    
    # Fallback: filename
    filename = os.path.splitext(os.path.basename(filepath))[0]
    match = state_pattern.search(filename.lower())
    if match:
        state_code = match.group(1).upper()
        mol_id = state_pattern.sub('', filename).strip('_-')
        return mol_id, state_code
    
    return filename, "unknown"


if __name__ == "__main__":
    main()
