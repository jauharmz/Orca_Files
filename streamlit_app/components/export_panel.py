"""
Export Panel Component - Complete Implementation

Features:
- Molecule + State selection for export
- CSV export (scalar data)
- Excel export (multiple sheets with state)
- JSON export (complete data)
- Download buttons
"""

import streamlit as st
import pandas as pd
import numpy as np
import json
from io import BytesIO
from typing import List


def render_export_panel(df: pd.DataFrame):
    """Render export panel with molecule+state selection."""
    
    st.subheader("📤 Export Data")
    
    if df.empty:
        st.warning("No data available to export")
        return
    
    # Build molecule options
    mol_options = []
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        if state and str(state) != "nan":
            label = f"{mol_id} [{state}]"
        else:
            label = mol_id
        mol_options.append({"label": label, "idx": idx})
    
    unique_labels = list(dict.fromkeys([m["label"] for m in mol_options]))
    
    # Selection - default all for export
    selected = st.multiselect(
        "Select Data to Export",
        unique_labels,
        default=unique_labels,
        key="export_mol_select"
    )
    
    if not selected:
        st.warning("Select data to export")
        return
    
    # Get selected data
    selected_indices = [m["idx"] for m in mol_options if m["label"] in selected]
    selected_df = df.loc[selected_indices]
    
    st.info(f"📊 {len(selected_df)} records selected for export")
    
    # Export tabs
    export_tabs = st.tabs(["📋 CSV", "📊 Excel", "📄 JSON", "🎨 HTML Report"])
    
    with export_tabs[0]:
        render_csv_export(selected_df)
    with export_tabs[1]:
        render_excel_export(selected_df)
    with export_tabs[2]:
        render_json_export(selected_df)
    with export_tabs[3]:
        render_html_export(selected_df)


def render_csv_export(df: pd.DataFrame):
    """Export scalar data to CSV."""
    
    st.markdown("##### CSV Export (Scalar Data)")
    
    # Select only scalar columns
    scalar_cols = [
        "molecule_id", "optimized_state", "filename", "smiles", "charge", "multiplicity",
        "gibbs_Eh", "single_point_Eh", "homo_energy", "lumo_energy", "homo_lumo_gap",
        "functional", "basis_set", "dispersion", "solvent",
        "is_optimization", "calc_class", "esd_type", "method_id"
    ]
    
    available_cols = [c for c in scalar_cols if c in df.columns]
    export_df = df[available_cols].copy()
    
    csv_data = export_df.to_csv(index=False)
    
    col1, col2 = st.columns([1, 3])
    with col1:
        st.download_button(
            "⬇️ Download CSV",
            csv_data,
            "orca_molecules.csv",
            "text/csv",
            key="csv_download"
        )
    with col2:
        st.caption(f"{len(export_df)} rows × {len(available_cols)} columns")
    
    with st.expander("Preview"):
        st.dataframe(export_df.head(10), use_container_width=True, hide_index=True)


def render_excel_export(df: pd.DataFrame):
    """Export to Excel with multiple sheets including state."""
    
    st.markdown("##### Excel Export (Multiple Sheets)")
    
    try:
        output = BytesIO()
        
        with pd.ExcelWriter(output, engine='openpyxl') as writer:
            # Sheet 1: Main molecule data
            scalar_cols = ["molecule_id", "optimized_state", "filename", "gibbs_Eh", 
                          "single_point_Eh", "homo_energy", "lumo_energy", "homo_lumo_gap",
                          "functional", "basis_set", "calc_class"]
            available = [c for c in scalar_cols if c in df.columns]
            df[available].to_excel(writer, sheet_name="Molecules", index=False)
            
            # Sheet 2: Cartesian coordinates
            all_coords = []
            for idx, row in df.iterrows():
                mol_id = row.get("molecule_id", "unknown")
                state = row.get("optimized_state", "")
                coords = row.get("cart_coords")
                if coords is not None and hasattr(coords, 'empty') and not coords.empty:
                    coords_copy = coords.copy()
                    coords_copy.insert(0, "molecule_id", mol_id)
                    coords_copy.insert(1, "optimized_state", state)
                    all_coords.append(coords_copy)
            if all_coords:
                pd.concat(all_coords, ignore_index=True).to_excel(writer, sheet_name="Coordinates", index=False)
            
            # Sheet 3: Orbitals
            all_orbitals = []
            for idx, row in df.iterrows():
                mol_id = row.get("molecule_id", "unknown")
                state = row.get("optimized_state", "")
                orbitals = row.get("orbitals")
                if orbitals is not None and hasattr(orbitals, 'empty') and not orbitals.empty:
                    orb_copy = orbitals.copy()
                    orb_copy.insert(0, "molecule_id", mol_id)
                    orb_copy.insert(1, "optimized_state", state)
                    all_orbitals.append(orb_copy)
            if all_orbitals:
                pd.concat(all_orbitals, ignore_index=True).to_excel(writer, sheet_name="Orbitals", index=False)
            
            # Sheet 4: IR Spectra
            all_ir = []
            for idx, row in df.iterrows():
                mol_id = row.get("molecule_id", "unknown")
                state = row.get("optimized_state", "")
                ir = row.get("ir")
                if ir is not None and hasattr(ir, 'empty') and not ir.empty:
                    ir_copy = ir.copy()
                    ir_copy.insert(0, "molecule_id", mol_id)
                    ir_copy.insert(1, "optimized_state", state)
                    all_ir.append(ir_copy)
            if all_ir:
                pd.concat(all_ir, ignore_index=True).to_excel(writer, sheet_name="IR_Spectra", index=False)
            
            # Sheet 5: TDDFT
            all_tddft = []
            for idx, row in df.iterrows():
                mol_id = row.get("molecule_id", "unknown")
                state = row.get("optimized_state", "")
                tddft = row.get("tddft_states")
                if tddft is not None and hasattr(tddft, 'empty') and not tddft.empty:
                    tddft_copy = tddft.copy()
                    tddft_copy.insert(0, "molecule_id", mol_id)
                    tddft_copy.insert(1, "optimized_state", state)
                    all_tddft.append(tddft_copy)
            if all_tddft:
                pd.concat(all_tddft, ignore_index=True).to_excel(writer, sheet_name="TDDFT", index=False)
        
        output.seek(0)
        
        st.download_button(
            "⬇️ Download Excel",
            output,
            "orca_data.xlsx",
            "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            key="excel_download"
        )
        
        st.caption("Sheets: Molecules, Coordinates, Orbitals, IR_Spectra, TDDFT")
        
    except ImportError:
        st.error("openpyxl not installed. Run: pip install openpyxl")


def render_json_export(df: pd.DataFrame):
    """Export complete data to JSON."""
    
    st.markdown("##### JSON Export (Complete Data)")
    
    def make_serializable(val):
        if isinstance(val, pd.DataFrame):
            if val.empty:
                return None
            return val.to_dict(orient='records')
        if pd.isna(val):
            return None
        if isinstance(val, (np.integer, np.floating)):
            return val.item()
        return val
    
    records = []
    for idx, row in df.iterrows():
        record = {}
        for col in df.columns:
            record[col] = make_serializable(row[col])
        records.append(record)
    
    json_str = json.dumps(records, indent=2, default=str)
    
    st.download_button(
        "⬇️ Download JSON",
        json_str,
        "orca_data.json",
        "application/json",
        key="json_download"
    )
    
    with st.expander("Preview (first record)"):
        if records:
            st.json(records[0])


def render_html_export(df: pd.DataFrame):
    """Export comprehensive interactive HTML report."""
    
    st.markdown("##### 🎨 Interactive HTML Report")
    st.info("📄 Generates a **self-contained HTML file** with interactive 3D molecules, spectra, and data tables.")
    
    # Initialize session state keys
    if 'html_report_data' not in st.session_state:
        st.session_state.html_report_data = None
    if 'html_report_name' not in st.session_state:
        st.session_state.html_report_name = "orca_report.html"
    if 'html_save_success' not in st.session_state:
        st.session_state.html_save_success = None
    if 'html_save_path' not in st.session_state:
        st.session_state.html_save_path = None
    if 'html_save_error' not in st.session_state:
        st.session_state.html_save_error = None
    
    # Export options
    with st.expander("⚙️ Export Options", expanded=True):
        c1, c2, c3 = st.columns(3)
        with c1:
            include_3d = st.checkbox("3D Molecule Viewer", True, key="html_3d")
            include_spectra = st.checkbox("Spectra Charts", True, key="html_spectra")
        with c2:
            include_energy = st.checkbox("Energy Diagram", True, key="html_energy")
            include_orbitals = st.checkbox("Orbital Charts", True, key="html_orbitals")
        with c3:
            dark_theme = st.checkbox("Dark Theme", False, key="html_dark")
            compact_mode = st.checkbox("Compact Mode", False, key="html_compact")
    
    # Report filename
    report_name = st.text_input("Report Filename", value="orca_report.html", key="html_filename_input")
    if not report_name.endswith(".html"):
        report_name += ".html"
    
    # Save path configuration
    from pathlib import Path
    import tempfile
    
    # Default save locations (in order of preference)
    default_paths = []
    
    # 1. Current working directory
    cwd = Path.cwd()
    default_paths.append(str(cwd))
    
    # 2. User home directory
    home = Path.home()
    default_paths.append(str(home))
    
    # 3. Downloads folder
    downloads = home / "Downloads"
    if downloads.exists():
        default_paths.append(str(downloads))
    
    # 4. Temp directory as fallback
    default_paths.append(tempfile.gettempdir())
    
    # Remove duplicates while preserving order
    seen = set()
    unique_paths = []
    for p in default_paths:
        if p not in seen:
            seen.add(p)
            unique_paths.append(p)
    
    save_directory = st.selectbox(
        "💾 Save Directory",
        options=unique_paths,
        index=0,
        key="html_save_dir",
        help="Select where to save the HTML file locally"
    )
    
    # Option to enter custom path
    custom_path = st.text_input(
        "Or enter custom path (optional)",
        value="",
        key="html_custom_path",
        help="Leave empty to use the selected directory above"
    )
    
    if custom_path.strip():
        save_directory = custom_path.strip()
    
    # Generate button
    col1, col2 = st.columns(2)
    with col1:
        generate_btn = st.button("🚀 Generate Report", type="primary", use_container_width=True)
    with col2:
        clear_btn = st.button("🗑️ Clear", use_container_width=True)
    
    if clear_btn:
        st.session_state.html_report_data = None
        st.session_state.html_report_name = "orca_report.html"
        st.session_state.html_save_success = None
        st.session_state.html_save_path = None
        st.session_state.html_save_error = None
        st.rerun()
    
    if generate_btn:
        try:
            with st.spinner("Generating comprehensive report (this may take a moment)..."):
                html = generate_html_report(
                    df, 
                    include_3d=include_3d,
                    include_spectra=include_spectra,
                    include_energy=include_energy,
                    include_orbitals=include_orbitals,
                    dark_theme=dark_theme,
                    compact_mode=compact_mode
                )
            
            # Store in session state
            st.session_state.html_report_data = html
            st.session_state.html_report_name = report_name
            
            # Try to save to disk
            save_path = Path(save_directory) / report_name
            try:
                # Ensure directory exists
                save_path.parent.mkdir(parents=True, exist_ok=True)
                with open(save_path, "w", encoding="utf-8") as f:
                    f.write(html)
                st.session_state.html_save_success = True
                st.session_state.html_save_path = str(save_path)
                st.session_state.html_save_error = None
            except PermissionError:
                st.session_state.html_save_success = False
                st.session_state.html_save_error = f"Permission denied: Cannot write to {save_path}"
            except Exception as save_err:
                st.session_state.html_save_success = False
                st.session_state.html_save_error = f"Failed to save: {save_err}"
            
            st.rerun()  # Rerun to show results
                
        except Exception as e:
            st.error(f"❌ Failed to generate report: {e}")
            import traceback
            st.code(traceback.format_exc())
    
    # Show download button if report exists
    if st.session_state.html_report_data:
        html = st.session_state.html_report_data
        name = st.session_state.html_report_name
        
        st.success(f"✅ Report generated! Size: {len(html.encode('utf-8')) / 1024:.1f} KB")
        
        # Show save status
        if st.session_state.html_save_success is True:
            st.caption(f"💾 Saved to: `{st.session_state.html_save_path}`")
        elif st.session_state.html_save_success is False:
            st.warning(f"⚠️ {st.session_state.html_save_error}")
            st.caption("Use the download button below instead.")
        
        # Always show download button - this is the primary download method
        st.download_button(
            label="📥 Download HTML Report",
            data=html,
            file_name=name,
            mime="text/html",
            type="primary",
            use_container_width=True,
            key="html_dl_btn"
        )



def generate_html_report(
    df: pd.DataFrame, 
    include_3d: bool = True,
    include_spectra: bool = True,
    include_energy: bool = True,
    include_orbitals: bool = True,
    dark_theme: bool = False,
    compact_mode: bool = False
) -> str:
    """Generate comprehensive research paper style HTML report."""
    
    # Theme colors
    if dark_theme:
        bg_primary = "#1a1a2e"
        bg_secondary = "#16213e"
        bg_card = "#1f3460"
        text_primary = "#eaeaea"
        text_secondary = "#a0a0b0"
        accent = "#00d4aa"
        accent2 = "#7c3aed"
    else:
        bg_primary = "#f5f7fa"
        bg_secondary = "#ffffff"
        bg_card = "#ffffff"
        text_primary = "#1a1a2e"
        text_secondary = "#64748b"
        accent = "#0ea5e9"
        accent2 = "#8b5cf6"
    
    # Calculate statistics for the report
    n_molecules = df["molecule_id"].nunique() if "molecule_id" in df.columns else 0
    n_records = len(df)
    n_states = df["optimized_state"].nunique() if "optimized_state" in df.columns else 0
    n_with_energy = df["gibbs_Eh"].notna().sum() if "gibbs_Eh" in df.columns else 0
    n_with_ir = df["ir"].apply(lambda x: x is not None and hasattr(x, 'empty') and not x.empty).sum() if "ir" in df.columns else 0
    n_with_raman = df["raman"].apply(lambda x: x is not None and hasattr(x, 'empty') and not x.empty).sum() if "raman" in df.columns else 0
    
    # Get unique methods
    methods = df["method_id"].unique().tolist() if "method_id" in df.columns else ["Unknown"]
    methods_str = ", ".join([str(m) for m in methods[:5]]) + ("..." if len(methods) > 5 else "")
    
    # Energy range
    energies = df["gibbs_Eh"].dropna() if "gibbs_Eh" in df.columns else pd.Series()
    if len(energies) > 0:
        e_min, e_max = energies.min(), energies.max()
        e_range_str = f"{e_min:.4f} to {e_max:.4f} Eh"
    else:
        e_range_str = "N/A"
    
    # HOMO-LUMO gap statistics
    gaps = df["homo_lumo_gap"].dropna() if "homo_lumo_gap" in df.columns else pd.Series()
    if len(gaps) > 0:
        gap_avg = gaps.mean()
        gap_std = gaps.std()
        gap_str = f"{gap_avg:.2f} ± {gap_std:.2f} eV"
    else:
        gap_str = "N/A"
    
    # Prepare molecule data as JSON
    molecules_json = prepare_molecules_json(df)
    
    # Current timestamp
    timestamp = pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")
    date_only = pd.Timestamp.now().strftime("%B %d, %Y")
    
    # Build HTML
    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ORCA Computational Chemistry Report</title>
    
    <!-- 3Dmol.js for 3D molecules -->
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    
    <!-- Plotly for charts -->
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    
    <!-- Google Fonts -->
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
    
    <style>
        :root {{
            --bg-primary: {bg_primary};
            --bg-secondary: {bg_secondary};
            --bg-card: {bg_card};
            --text-primary: {text_primary};
            --text-secondary: {text_secondary};
            --accent: {accent};
            --accent2: {accent2};
        }}
        
        * {{ box-sizing: border-box; margin: 0; padding: 0; }}
        
        body {{
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
            background: var(--bg-primary);
            color: var(--text-primary);
            line-height: 1.7;
            font-size: 15px;
        }}
        
        .container {{ max-width: 1200px; margin: 0 auto; padding: 40px 20px; }}
        
        /* Header */
        .report-header {{
            text-align: center;
            padding: 60px 40px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 50%, #f472b6 100%);
            color: white;
            border-radius: 20px;
            margin-bottom: 40px;
            box-shadow: 0 20px 60px rgba(102, 126, 234, 0.3);
        }}
        
        .report-header h1 {{
            font-size: 2.5em;
            font-weight: 700;
            margin-bottom: 10px;
            letter-spacing: -0.5px;
        }}
        
        .report-header .subtitle {{
            font-size: 1.2em;
            opacity: 0.9;
            margin-bottom: 20px;
        }}
        
        .report-header .meta {{
            font-size: 0.9em;
            opacity: 0.8;
        }}
        
        /* Navigation */
        .toc {{
            background: var(--bg-card);
            border-radius: 16px;
            padding: 30px;
            margin-bottom: 40px;
            box-shadow: 0 4px 20px rgba(0,0,0,0.08);
        }}
        
        .toc h2 {{ 
            color: var(--accent);
            margin-bottom: 20px;
            font-size: 1.3em;
        }}
        
        .toc-list {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 10px;
            list-style: none;
        }}
        
        .toc-list a {{
            color: var(--text-primary);
            text-decoration: none;
            padding: 12px 16px;
            border-radius: 8px;
            display: block;
            transition: all 0.2s;
            border-left: 3px solid transparent;
        }}
        
        .toc-list a:hover {{
            background: var(--bg-primary);
            border-left-color: var(--accent);
            transform: translateX(5px);
        }}
        
        /* Sections */
        .section {{
            background: var(--bg-card);
            border-radius: 16px;
            padding: 40px;
            margin-bottom: 30px;
            box-shadow: 0 4px 20px rgba(0,0,0,0.05);
        }}
        
        .section h2 {{
            color: var(--accent);
            font-size: 1.6em;
            margin-bottom: 25px;
            padding-bottom: 15px;
            border-bottom: 2px solid var(--accent);
            display: flex;
            align-items: center;
            gap: 12px;
        }}
        
        .section h3 {{
            color: var(--text-primary);
            font-size: 1.2em;
            margin: 25px 0 15px 0;
        }}
        
        .section p {{
            color: var(--text-secondary);
            margin-bottom: 15px;
        }}
        
        /* Executive Summary Box */
        .summary-box {{
            background: linear-gradient(135deg, rgba(14, 165, 233, 0.1) 0%, rgba(139, 92, 246, 0.1) 100%);
            border-left: 4px solid var(--accent);
            padding: 25px;
            border-radius: 0 12px 12px 0;
            margin: 20px 0;
        }}
        
        .summary-box p {{
            color: var(--text-primary);
            font-size: 1.05em;
        }}
        
        /* Stats Grid */
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 30px 0;
        }}
        
        .stat-card {{
            background: var(--bg-primary);
            padding: 25px;
            border-radius: 12px;
            text-align: center;
            border: 1px solid rgba(0,0,0,0.05);
            transition: transform 0.2s, box-shadow 0.2s;
        }}
        
        .stat-card:hover {{
            transform: translateY(-3px);
            box-shadow: 0 10px 30px rgba(0,0,0,0.1);
        }}
        
        .stat-value {{
            font-size: 2.2em;
            font-weight: 700;
            color: var(--accent);
            font-family: 'JetBrains Mono', monospace;
        }}
        
        .stat-label {{
            color: var(--text-secondary);
            font-size: 0.9em;
            margin-top: 5px;
            font-weight: 500;
        }}
        
        /* Methods Table */
        .methods-table {{
            width: 100%;
            margin: 20px 0;
            border-collapse: collapse;
        }}
        
        .methods-table th, .methods-table td {{
            padding: 15px;
            text-align: left;
            border-bottom: 1px solid rgba(0,0,0,0.1);
        }}
        
        .methods-table th {{
            background: var(--bg-primary);
            font-weight: 600;
            color: var(--text-primary);
        }}
        
        .methods-table tr:hover td {{
            background: var(--bg-primary);
        }}
        
        /* 3D Viewer */
        .viewer-container {{
            background: {"#0f0f1a" if dark_theme else "#f0f4f8"};
            border-radius: 12px;
            padding: 20px;
            margin: 20px 0;
            position: relative;
            z-index: 1;
        }}
        
        #viewer-3d {{
            width: 100%;
            height: 500px;
            border-radius: 8px;
            border: 1px solid rgba(0,0,0,0.1);
            position: relative;
            z-index: 1;
        }}
        
        .controls {{
            display: flex;
            gap: 10px;
            flex-wrap: wrap;
            margin-bottom: 20px;
            align-items: center;
        }}
        
        select, button {{
            padding: 10px 18px;
            border-radius: 8px;
            border: 1px solid rgba(0,0,0,0.1);
            background: var(--bg-card);
            color: var(--text-primary);
            font-size: 14px;
            font-weight: 500;
            cursor: pointer;
            transition: all 0.2s;
        }}
        
        button:hover {{ 
            background: var(--accent); 
            color: white;
            transform: translateY(-2px);
        }}
        
        button.active {{ 
            background: var(--accent); 
            color: white;
        }}
        
        /* Chart Container */
        .chart-container {{ 
            width: 100%; 
            height: 450px;
            margin: 20px 0;
        }}
        
        /* Tabs */
        .tabs {{
            display: flex;
            gap: 5px;
            margin-bottom: 0;
            border-bottom: 2px solid var(--bg-primary);
        }}
        
        .tab {{
            padding: 12px 24px;
            border-radius: 8px 8px 0 0;
            cursor: pointer;
            font-weight: 500;
            transition: all 0.2s;
            background: var(--bg-primary);
        }}
        
        .tab:hover {{ background: var(--accent); color: white; }}
        .tab.active {{ background: var(--accent); color: white; }}
        
        .tab-content {{ display: none; padding-top: 20px; }}
        .tab-content.active {{ display: block; }}
        
        /* Data Table */
        .data-table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 14px;
            margin-top: 20px;
        }}
        
        .data-table th, .data-table td {{
            padding: 12px 15px;
            text-align: left;
            border-bottom: 1px solid rgba(0,0,0,0.08);
        }}
        
        .data-table th {{
            background: var(--bg-primary);
            font-weight: 600;
            color: var(--text-primary);
            position: sticky;
            top: 0;
        }}
        
        .data-table tbody tr:hover {{
            background: var(--bg-primary);
        }}
        
        .data-table code {{
            font-family: 'JetBrains Mono', monospace;
            font-size: 0.9em;
            background: var(--bg-primary);
            padding: 2px 6px;
            border-radius: 4px;
        }}
        
        /* Molecule Info Cards */
        .mol-info-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(140px, 1fr));
            gap: 15px;
            margin-top: 20px;
        }}
        
        /* Interactive Settings Panel */
        .settings-panel {{
            background: var(--bg-primary);
            border-radius: 10px;
            padding: 15px;
            margin: 15px 0;
            border: 1px solid rgba(0,0,0,0.1);
        }}
        
        .settings-panel summary {{
            cursor: pointer;
            font-weight: 600;
            color: var(--accent);
            padding: 5px 0;
            list-style: none;
        }}
        
        .settings-panel summary::-webkit-details-marker {{ display: none; }}
        .settings-panel summary::before {{ content: "⚙️ "; }}
        
        .settings-row {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 15px;
            margin-top: 15px;
        }}
        
        .setting-group {{
            display: flex;
            flex-direction: column;
            gap: 5px;
        }}
        
        .setting-group label {{
            font-size: 0.85em;
            font-weight: 500;
            color: var(--text-secondary);
        }}
        
        .setting-group input[type="range"] {{
            width: 100%;
            accent-color: var(--accent);
        }}
        
        .setting-group input[type="checkbox"] {{
            width: 18px;
            height: 18px;
            accent-color: var(--accent);
        }}
        
        .setting-group select {{
            padding: 8px 12px;
            border-radius: 6px;
            border: 1px solid rgba(0,0,0,0.1);
            background: var(--bg-card);
            color: var(--text-primary);
            font-size: 14px;
        }}
        
        .range-value {{
            font-family: 'JetBrains Mono', monospace;
            font-size: 0.9em;
            color: var(--accent);
            font-weight: 600;
        }}
        
        /* Multi-Select Dropdown */
        .multiselect-container {{
            position: relative;
            width: 100%;
        }}
        
        .multiselect-dropdown {{
            width: 100%;
            padding: 10px 14px;
            border-radius: 8px;
            border: 1px solid rgba(0,0,0,0.15);
            background: var(--bg-card);
            color: var(--text-primary);
            font-size: 14px;
            cursor: pointer;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }}
        
        .multiselect-dropdown:hover {{
            border-color: var(--accent);
        }}
        
        .multiselect-options {{
            position: absolute;
            top: 100%;
            left: 0;
            right: 0;
            max-height: 250px;
            overflow-y: auto;
            background: var(--bg-card);
            border: 1px solid rgba(0,0,0,0.15);
            border-radius: 8px;
            box-shadow: 0 8px 24px rgba(0,0,0,0.15);
            z-index: 100;
            display: none;
            margin-top: 4px;
        }}
        
        .multiselect-options.show {{
            display: block;
        }}
        
        .multiselect-option {{
            padding: 10px 14px;
            cursor: pointer;
            display: flex;
            align-items: center;
            gap: 10px;
            transition: background 0.15s;
        }}
        
        .multiselect-option:hover {{
            background: var(--bg-primary);
        }}
        
        .multiselect-option.selected {{
            background: rgba(14, 165, 233, 0.1);
        }}
        
        .multiselect-option input[type="checkbox"] {{
            width: 16px;
            height: 16px;
            accent-color: var(--accent);
        }}
        
        .multiselect-pills {{
            display: flex;
            flex-wrap: wrap;
            gap: 6px;
            margin-top: 10px;
        }}
        
        .pill {{
            display: inline-flex;
            align-items: center;
            gap: 6px;
            padding: 4px 10px;
            background: var(--accent);
            color: white;
            border-radius: 20px;
            font-size: 12px;
            font-weight: 500;
        }}
        
        .pill .remove {{
            cursor: pointer;
            opacity: 0.8;
            font-size: 14px;
        }}
        
        .pill .remove:hover {{
            opacity: 1;
        }}
        
        .multiselect-actions {{
            display: flex;
            gap: 8px;
            padding: 8px 14px;
            border-bottom: 1px solid rgba(0,0,0,0.1);
        }}
        
        .multiselect-actions button {{
            padding: 4px 10px;
            font-size: 12px;
            border-radius: 4px;
        }}
        
        /* Multi-Select Grid Layout */
        .multiselect-grid {{
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
            margin-top: 10px;
            max-height: 200px;
            overflow-y: auto;
            padding: 10px;
            background: var(--bg-primary);
            border-radius: 8px;
            border: 1px solid rgba(0,0,0,0.1);
        }}
        
        .multiselect-item {{
            display: inline-flex;
            align-items: center;
            gap: 6px;
            padding: 6px 12px;
            background: var(--bg-card);
            border: 1px solid rgba(0,0,0,0.08);
            border-radius: 20px;
            font-size: 13px;
            cursor: pointer;
            transition: all 0.15s;
        }}
        
        .multiselect-item:hover {{
            background: var(--accent);
            color: white;
            border-color: var(--accent);
        }}
        
        .multiselect-item input {{
            width: 14px;
            height: 14px;
            accent-color: var(--accent);
        }}
        
        .multiselect-item:has(input:checked) {{
            background: var(--accent);
            color: white;
            border-color: var(--accent);
        }}
        
        .multiselect-actions-bar {{
            display: flex;
            gap: 8px;
            margin-top: 10px;
        }}
        
        .multiselect-actions-bar button {{
            padding: 4px 10px;
            font-size: 11px;
            border-radius: 4px;
            background: var(--bg-card);
        }}
        
        /* Properties Tabs (Sub-tabs) */
        .prop-tabs {{
            display: flex;
            gap: 4px;
            margin-bottom: 15px;
            border-bottom: 1px solid rgba(0,0,0,0.1);
            padding-bottom: 0;
        }}
        
        .prop-tab {{
            padding: 8px 16px;
            border-radius: 6px 6px 0 0;
            cursor: pointer;
            font-size: 13px;
            font-weight: 500;
            background: transparent;
            color: var(--text-secondary);
            border: none;
            transition: all 0.2s;
        }}
        
        .prop-tab:hover {{
            background: var(--bg-primary);
            color: var(--text-primary);
        }}
        
        .prop-tab.active {{
            background: var(--accent);
            color: white;
        }}
        
        .prop-content {{
            display: none;
        }}
        
        .prop-content.active {{
            display: block;
        }}
        
        /* Coordinate/Data Tables */
        .data-scroll {{
            max-height: 400px;
            overflow-y: auto;
            border-radius: 8px;
            margin-top: 15px;
        }}
        
        /* Compact Data Table */
        .compact-table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 13px;
        }}
        
        .compact-table th, .compact-table td {{
            padding: 8px 12px;
            text-align: left;
            border-bottom: 1px solid rgba(0,0,0,0.08);
        }}
        
        .compact-table th {{
            background: var(--bg-primary);
            font-weight: 600;
            position: sticky;
            top: 0;
        }}
        
        .compact-table tr:hover td {{
            background: var(--bg-primary);
        }}
        
        /* Range Input with Dual Handles */
        .dual-range {{
            display: flex;
            align-items: center;
            gap: 10px;
        }}
        
        .dual-range input[type="number"] {{
            width: 70px;
            padding: 6px 8px;
            border-radius: 4px;
            border: 1px solid rgba(0,0,0,0.15);
            background: var(--bg-card);
            color: var(--text-primary);
            font-size: 13px;
            text-align: center;
        }}
        
        /* Inline Checkbox Label */
        .checkbox-inline {{
            display: flex;
            align-items: center;
            gap: 8px;
            cursor: pointer;
        }}
        
        .checkbox-inline input {{
            width: 16px;
            height: 16px;
        }}
        
        /* Footer */
        footer {{
            text-align: center;
            padding: 40px 20px;
            color: var(--text-secondary);
            font-size: 13px;
            border-top: 1px solid rgba(0,0,0,0.1);
            margin-top: 40px;
        }}
        
        footer a {{ color: var(--accent); text-decoration: none; }}
        
        /* Print Styles */
        @media print {{
            .controls, .tabs {{ display: none !important; }}
            .section {{ break-inside: avoid; }}
        }}
        
        /* Responsive */
        @media (max-width: 768px) {{
            .container {{ padding: 20px 10px; }}
            .section {{ padding: 20px 15px; }}
            .report-header {{ padding: 30px 15px; }}
            .report-header h1 {{ font-size: 1.5em; }}
            .stats-grid {{ grid-template-columns: repeat(2, 1fr); gap: 12px; }}
            .mol-info-grid {{ grid-template-columns: repeat(2, 1fr); gap: 8px; }}
            .mol-info-grid .stat-card {{ padding: 15px 10px; }}
            .mol-info-grid .stat-value {{ font-size: 1.4em; }}
            .mol-info-grid .stat-label {{ font-size: 0.8em; }}
            .controls {{ flex-direction: column; gap: 8px; align-items: stretch; }}
            .controls button, .controls select {{ width: 100%; justify-content: center; }}
            #viewer-3d {{ height: 350px; }}
            .chart-container {{ height: 350px; }}
            .tabs {{ flex-wrap: wrap; gap: 4px; }}
            .tab {{ flex: 1 1 auto; text-align: center; padding: 10px 8px; font-size: 13px; min-width: 80px; }}
            table {{ font-size: 12px; }}
            td, th {{ padding: 8px 6px; }}
            .settings-row {{ grid-template-columns: 1fr 1fr !important; gap: 12px; }}
            .setting-group {{ margin-bottom: 10px; }}
            .prop-tabs {{ flex-wrap: wrap; gap: 4px; }}
            .prop-tab {{ font-size: 11px; padding: 6px 10px; flex: 1 1 auto; text-align: center; }}
            .multiselect-pills {{ gap: 4px; }}
            .pill {{ font-size: 11px; padding: 3px 8px; }}
            .table-scroll {{ overflow-x: auto; -webkit-overflow-scrolling: touch; }}
            .viewer-container {{ padding: 12px; }}
        }}
        
        @media (max-width: 480px) {{
            .report-header h1 {{ font-size: 1.2em; }}
            .report-header .subtitle {{ font-size: 1em; }}
            .section {{ padding: 12px 8px; }}
            .section h2 {{ font-size: 1.3em; gap: 8px; }}
            #viewer-3d {{ height: 260px; }}
            .chart-container {{ height: 260px; }}
            .tab {{ padding: 8px 4px; font-size: 10px; min-width: 55px; }}
            .mol-info-grid {{ grid-template-columns: repeat(2, 1fr); gap: 6px; }}
            .mol-info-grid .stat-card {{ padding: 10px 6px; }}
            .mol-info-grid .stat-value {{ font-size: 1.1em; }}
            .mol-info-grid .stat-label {{ font-size: 0.7em; }}
            .stats-grid {{ gap: 8px; }}
            .stats-grid .stat-card {{ padding: 12px 8px; }}
            .stats-grid .stat-value {{ font-size: 1.6em; }}
            .stats-grid .stat-label {{ font-size: 0.8em; }}
            .settings-row {{ grid-template-columns: 1fr !important; }}
            .controls button {{ padding: 8px 12px; font-size: 13px; }}
            .prop-tabs {{ gap: 2px; }}
            .prop-tab {{ font-size: 10px; padding: 5px 8px; }}
            .dual-range {{ flex-wrap: wrap; }}
            .dual-range input[type="number"] {{ width: 60px; }}
        }}
        
        /* Table scroll wrapper */
        .table-scroll {{
            overflow-x: auto;
            -webkit-overflow-scrolling: touch;
            margin: 0 -10px;
            padding: 0 10px;
        }}
        
        /* Table styles to prevent truncation */
        table td, table th {{
            word-break: break-word;
            overflow-wrap: break-word;
            max-width: 300px;
        }}
        
        table td.smiles-cell {{
            font-family: 'JetBrains Mono', monospace;
            font-size: 11px;
            max-width: 200px;
            word-break: break-all;
        }}
    </style>
</head>
<body>
    <div class="container">
        <!-- Report Header -->
        <header class="report-header">
            <h1>🔬 Computational Chemistry Report</h1>
            <p class="subtitle">ORCA Quantum Chemical Calculations Analysis</p>
            <p class="meta">Generated: {date_only} | Molecules: {n_molecules} | Dataset: {n_records} calculations</p>
        </header>
        
        <!-- Table of Contents -->
        <nav class="toc">
            <h2>📑 Contents</h2>
            <ul class="toc-list">
                <li><a href="#executive-summary">1. Executive Summary</a></li>
                <li><a href="#methodology">2. Computational Methods</a></li>
                <li><a href="#molecular-structures">3. Molecular Structures</a></li>
                <li><a href="#electronic-structure">4. Electronic Structure</a></li>
                <li><a href="#energy-analysis">5. Energy Analysis</a></li>
                <li><a href="#spectra-analysis">6. Vibrational Analysis</a></li>
                <li><a href="#data-appendix">7. Complete Data Appendix</a></li>
            </ul>
        </nav>
        
        <!-- 1. Executive Summary -->
        <section id="executive-summary" class="section">
            <h2>📋 1. Executive Summary</h2>
            
            <div class="summary-box">
                <p>
                    This report presents comprehensive computational chemistry results for <strong>{n_molecules} molecular system{"s" if n_molecules != 1 else ""}</strong> 
                    analyzed using ORCA quantum chemistry software. A total of <strong>{n_records} calculations</strong> were performed 
                    covering <strong>{n_states} electronic state{"s" if n_states != 1 else ""}</strong> (ground and excited states).
                </p>
            </div>
            
            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-value">{n_molecules}</div>
                    <div class="stat-label">Molecules Studied</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{n_records}</div>
                    <div class="stat-label">Total Calculations</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{n_states}</div>
                    <div class="stat-label">Electronic States</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{n_with_ir}</div>
                    <div class="stat-label">IR Spectra</div>
                </div>
            </div>
            
            <h3>Key Findings</h3>
            <ul style="color: var(--text-secondary); margin-left: 20px;">
                <li><strong>Energy Range:</strong> {e_range_str}</li>
                <li><strong>Average HOMO-LUMO Gap:</strong> {gap_str}</li>
                <li><strong>Computational Methods:</strong> {methods_str}</li>
                <li><strong>Vibrational Data:</strong> {n_with_ir} IR spectra, {n_with_raman} Raman spectra available</li>
            </ul>
        </section>
        
        <!-- 2. Methodology -->
        <section id="methodology" class="section">
            <h2>⚗️ 2. Computational Methods</h2>
            
            <p>
                All calculations were performed using the <strong>ORCA</strong> quantum chemistry program package. 
                The following computational methods and parameters were employed:
            </p>
            
            <table class="methods-table">
                <thead>
                    <tr>
                        <th>Parameter</th>
                        <th>Value</th>
                        <th>Description</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td><strong>Electronic Structure Method</strong></td>
                        <td><code>{methods[0] if methods else "DFT"}</code></td>
                        <td>Primary computational method</td>
                    </tr>
                    <tr>
                        <td><strong>Number of Calculations</strong></td>
                        <td>{n_records}</td>
                        <td>Total geometry optimizations and single points</td>
                    </tr>
                    <tr>
                        <td><strong>Electronic States</strong></td>
                        <td>{", ".join([str(s) for s in df["optimized_state"].unique().tolist()[:5]]) if "optimized_state" in df.columns else "S0"}</td>
                        <td>Ground and excited state calculations</td>
                    </tr>
                    <tr>
                        <td><strong>Vibrational Analysis</strong></td>
                        <td>{"Frequency calculations performed" if n_with_ir > 0 else "Not included"}</td>
                        <td>IR and Raman spectroscopy prediction</td>
                    </tr>
                </tbody>
            </table>
        </section>
        
        <!-- 3. Molecular Structures -->
        <section id="molecular-structures" class="section">
            <h2>🧬 3. Molecular Structures</h2>
            
            <p>
                Interactive 3D visualization of the optimized molecular geometries. Select a molecule from the dropdown 
                and use the style buttons to change the representation.
            </p>
            
            <div class="viewer-container">
                <div class="controls" style="margin-bottom: 10px;">
                    <label><strong>Select Molecule:</strong></label>
                    <select id="molecule-select" onchange="switchMolecule()">
'''
    
    # Add molecule options
    for i, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        smiles = row.get("smiles", "")
        label = f"{mol_id} [{state}]" if state and str(state) != "nan" else mol_id
        html += f'                        <option value="{i}">{label}</option>\n'
    
    html += f'''                </select>
                </div>
                
                <div class="tabs">
                    <div class="tab active" onclick="showMolTab('3d')">🔮 3D View</div>
                    <div class="tab" onclick="showMolTab('2d')">📐 2D View</div>
                    <div class="tab" onclick="showMolTab('props')">📊 Properties</div>
                </div>

                <!-- 3D View Content -->
                <div id="mol-view-3d" class="tab-content active">
                    <div class="controls">
'''
    
    if include_3d:
        html += '''                    <button onclick="setStyle('stick')" id="btn-stick">Stick</button>
                    <button onclick="setStyle('sphere')" id="btn-sphere">Spacefill</button>
                    <button onclick="setStyle('ballstick')" id="btn-ballstick" class="active">Ball & Stick</button>
                    <button onclick="toggleSpin()" id="btn-spin">🔄 Spin</button>
                    <button onclick="toggleLabels()" id="btn-labels">🏷️ Labels</button>
'''
    
    html += '''                    </div>
'''
    
    if include_3d:
        html += f'''
                    <!-- Settings Panel -->
                    <details class="settings-panel">
                        <summary>3D Visualization Settings</summary>
                        <div class="settings-row">
                            <div class="setting-group">
                                <label>Bond Radius</label>
                                <input type="range" id="bond-radius" min="0.05" max="0.30" step="0.01" value="0.12" onchange="updateBondRadius(this.value)">
                                <span class="range-value" id="bond-radius-val">0.12</span>
                            </div>
                            <div class="setting-group">
                                <label>Atom Scale</label>
                                <input type="range" id="atom-scale" min="0.1" max="1.0" step="0.05" value="0.25" onchange="updateAtomScale(this.value)">
                                <span class="range-value" id="atom-scale-val">0.25</span>
                            </div>
                            <div class="setting-group">
                                <label>Color Scheme</label>
                                <select id="color-scheme" onchange="updateColorScheme(this.value)">
                                    <option value="Jmol">Standard (Jmol)</option>
                                    <option value="rasmol">RasMol</option>
                                    <option value="default">Default</option>
                                </select>
                            </div>
                            <div class="setting-group">
                                <label>Background</label>
                                <select id="bg-color" onchange="updateBackground(this.value)">
                                    <option value="{"#1e1e2e" if dark_theme else "#f0f0f5"}">Default</option>
                                    <option value="#000000">Black</option>
                                    <option value="#ffffff">White</option>
                                    <option value="#1a1a2e">Dark Blue</option>
                                </select>
                            </div>
                        </div>
                    </details>
                    
                    <div id="viewer-3d" style="background: {"#1e1e2e" if dark_theme else "#f0f0f5"};"></div>
'''

    html += f'''
                </div>
                
                <!-- 2D View Content -->
                <div id="mol-view-2d" class="tab-content">
                    <!-- 2D Settings Panel -->
                    <details class="settings-panel">
                        <summary>2D Visualization Settings</summary>
                        <div class="settings-row">
                            <div class="setting-group">
                                <label>Image Size</label>
                                <input type="range" id="mol2d-size" min="300" max="700" step="50" value="450" onchange="update2DSettings()">
                                <span class="range-value" id="mol2d-size-val">450px</span>
                            </div>
                            <div class="setting-group">
                                <label>Background</label>
                                <select id="mol2d-bg" onchange="update2DSettings()">
                                    <option value="#ffffff">White</option>
                                    <option value="#f5f5f5">Light Gray</option>
                                    <option value="#000000">Black</option>
                                    <option value="#1a1a24">Dark</option>
                                </select>
                            </div>
                            <div class="setting-group">
                                <label>Color Palette</label>
                                <select id="mol2d-palette" onchange="update2DSettings()">
                                    <option value="standard">Standard (RDKit)</option>
                                    <option value="bw">Black & White</option>
                                    <option value="dark">Dark Mode</option>
                                </select>
                            </div>
                            <div class="setting-group">
                                <label>Bond Width</label>
                                <input type="range" id="mol2d-bondwidth" min="1.0" max="4.0" step="0.5" value="2.0" onchange="update2DSettings()">
                                <span class="range-value" id="mol2d-bondwidth-val">2.0</span>
                            </div>
                        </div>
                        <div class="settings-row">
                            <div class="setting-group">
                                <label class="checkbox-inline">
                                    <input type="checkbox" id="mol2d-indices" onchange="update2DSettings()">
                                    Show Atom Indices
                                </label>
                            </div>
                            <div class="setting-group">
                                <label class="checkbox-inline">
                                    <input type="checkbox" id="mol2d-stereo" checked onchange="update2DSettings()">
                                    Show Stereo Labels
                                </label>
                            </div>
                            <div class="setting-group">
                                <label class="checkbox-inline">
                                    <input type="checkbox" id="mol2d-kekulize" onchange="update2DSettings()">
                                    Kekulize Bonds
                                </label>
                            </div>
                            <div class="setting-group">
                                <label>Padding</label>
                                <input type="range" id="mol2d-padding" min="0" max="0.2" step="0.02" value="0.05" onchange="update2DSettings()">
                                <span class="range-value" id="mol2d-padding-val">0.05</span>
                            </div>
                        </div>
                        <div class="settings-row">
                            <div class="setting-group">
                                <label>Atom Font Size</label>
                                <input type="range" id="mol2d-fontsize" min="10" max="24" step="2" value="14" onchange="update2DSettings()">
                                <span class="range-value" id="mol2d-fontsize-val">14</span>
                            </div>
                        </div>
                    </details>
                    
                    <div id="mol2d-container" style="display: flex; justify-content: center; align-items: center; background: #ffffff; border-radius: 8px; padding: 20px; min-height: 400px;">
                        <img id="img-2d" src="" alt="2D Structure" style="max-width: 100%; max-height: 400px; display: none;">
                        <div id="no-2d-msg" style="color: #666;">No 2D image available (requires SMILES)</div>
                    </div>
                </div>
                
                <!-- Properties Tab Content -->
                <div id="mol-view-props" class="tab-content">
                    <div class="prop-tabs">
                        <div class="prop-tab active" onclick="showPropSubTab('coords')">📍 Coordinates</div>
                        <div class="prop-tab" onclick="showPropSubTab('properties')">⚡ Properties</div>
                        <div class="prop-tab" onclick="showPropSubTab('mulliken')">💫 Mulliken</div>
                        <div class="prop-tab" onclick="showPropSubTab('method')">🔧 Method</div>
                    </div>
                    
                    <!-- Coordinates Sub-Tab -->
                    <div id="prop-coords" class="prop-content active">
                        <div class="data-scroll">
                            <table class="compact-table" id="coords-table">
                                <thead>
                                    <tr><th>#</th><th>Atom</th><th>X (Å)</th><th>Y (Å)</th><th>Z (Å)</th></tr>
                                </thead>
                                <tbody id="coords-tbody"></tbody>
                            </table>
                        </div>
                    </div>
                    
                    <!-- Properties Sub-Tab -->
                    <div id="prop-properties" class="prop-content">
                        <div class="data-scroll">
                            <table class="compact-table" id="props-table">
                                <thead>
                                    <tr><th>Property</th><th>Value</th></tr>
                                </thead>
                                <tbody id="props-tbody"></tbody>
                            </table>
                        </div>
                    </div>
                    
                    <!-- Mulliken Sub-Tab -->
                    <div id="prop-mulliken" class="prop-content">
                        <div class="data-scroll">
                            <table class="compact-table" id="mulliken-table">
                                <thead>
                                    <tr><th>#</th><th>Atom</th><th>Charge</th><th>Spin</th></tr>
                                </thead>
                                <tbody id="mulliken-tbody"></tbody>
                            </table>
                        </div>
                    </div>
                    
                    <!-- Method Sub-Tab -->
                    <div id="prop-method" class="prop-content">
                        <div class="data-scroll">
                            <table class="compact-table" id="method-table">
                                <thead>
                                    <tr><th>Parameter</th><th>Value</th></tr>
                                </thead>
                                <tbody id="method-tbody"></tbody>
                            </table>
                        </div>
                    </div>
                </div>
            
                <!-- Molecule Info -->
                <div id="mol-info" class="mol-info-grid"></div>
            </div>
        </section>
'''
    
    # Note: Energy Analysis section is now only included once (after Vibrational Analysis) with full controls
    
    if include_orbitals:
        html += '''        
        <!-- 4. Electronic Structure -->
        <section id="electronic-structure" class="section">
            <h2>🔮 4. Electronic Structure</h2>
            
            <p>
                HOMO (Highest Occupied Molecular Orbital) and LUMO (Lowest Unoccupied Molecular Orbital) energy levels 
                provide insight into the electronic properties, reactivity, and optical characteristics of the molecules.
            </p>
            
            <!-- Orbital Settings -->
            <details class="settings-panel" open>
                <summary>Orbital Visualization Settings</summary>
                <div class="settings-row">
                    <div class="setting-group">
                        <label>Select Molecule</label>
                        <select id="orbital-mol-select" onchange="renderOrbitals(parseInt(this.value))">
                            <!-- Options injected by JS -->
                        </select>
                    </div>
                    <div class="setting-group">
                        <label>Orbitals to Show (+/- HOMO)</label>
                        <input type="range" id="orb-range-slider" min="3" max="20" value="10" onchange="updateOrbitalRange(this.value)">
                        <span class="range-value" id="orb-range-val">10</span>
                    </div>
                    <div class="setting-group">
                        <label>Gap Mode</label>
                        <select id="orb-gap-mode" onchange="renderOrbitals(currentMolIndex)">
                            <option value="HL">HOMO-LUMO</option>
                            <option value="SL">SOMO-LUMO</option>
                            <option value="SS">SOMO-SUMO</option>
                        </select>
                    </div>
                    <div class="setting-group">
                        <label class="checkbox-inline">
                            <input type="checkbox" id="orb-connectors" onchange="renderOrbitals(currentMolIndex)">
                            Connector Lines
                        </label>
                    </div>
                    <div class="setting-group">
                        <label class="checkbox-inline">
                            <input type="checkbox" id="orb-labels" checked onchange="renderOrbitals(currentMolIndex)">
                            Orbital Labels
                        </label>
                    </div>
                </div>
                <div class="setting-group" style="margin-top: 10px; width: 100%;">
                    <label><strong>Select Molecules for Comparison:</strong></label>
                    <div id="orbital-multiselect" class="multiselect-grid">
                        <!-- Checkboxes injected by JS -->
                    </div>
                </div>
            </details>
            <div id="orbital-chart" class="chart-container"></div>
            <p style="font-size: 12px; color: var(--text-secondary); margin-top: 10px;">
                <strong>Legend:</strong> Solid = Occupied, Dashed = Virtual, Dot-Dash = SOMO/SUMO. Bold lines = HOMO/LUMO.
            </p>
        </section>
'''

    
    if include_spectra:
        html += '''        
        <!-- 6. Vibrational Analysis -->
        <section id="spectra-analysis" class="section">
            <h2>📈 6. Vibrational Analysis</h2>
            
            <p>
                Infrared (IR) and Raman spectra provide fingerprints of molecular vibrations. 
                Use the settings panel to customize the visualization.
            </p>
            
            <!-- Spectra Settings -->
            <details class="settings-panel" open>
                <summary>Spectra Visualization Settings</summary>
                <div class="settings-row">
                    <div class="setting-group">
                        <label>Frequency Range (cm⁻¹)</label>
                        <div class="dual-range">
                            <input type="number" id="freq-min" value="400" min="0" max="4000" onchange="renderSpectra()">
                            <span>to</span>
                            <input type="number" id="freq-max" value="4000" min="0" max="5000" onchange="renderSpectra()">
                        </div>
                    </div>
                    <div class="setting-group">
                        <label>FWHM Broadening (cm⁻¹)</label>
                        <input type="range" id="fwhm-slider" min="1" max="100" value="20" onchange="updateSpectraBroadening(this.value)">
                        <span class="range-value" id="fwhm-val">20 cm⁻¹</span>
                    </div>
                    <div class="setting-group">
                        <label>Display Mode</label>
                        <select id="spectra-display-mode" onchange="renderSpectra()">
                            <option value="overlay">Overlay</option>
                            <option value="stacked">Stacked</option>
                        </select>
                    </div>
                    <div class="setting-group">
                        <label>Visualization Mode</label>
                        <select id="spectra-mode" onchange="updateSpectraMode(this.value)">
                             <option value="line">Line (Smoothed)</option>
                             <option value="stick">Stick (Discrete)</option>
                        </select>
                    </div>
                </div>
                <div class="settings-row">
                    <div class="setting-group">
                        <label class="checkbox-inline">
                            <input type="checkbox" id="normalize-spectra" checked onchange="renderSpectra()">
                            Normalize Intensity
                        </label>
                    </div>
                    <div class="setting-group">
                        <label class="checkbox-inline">
                            <input type="checkbox" id="show-regions" onchange="renderSpectra()">
                            Show Region Boundaries
                        </label>
                    </div>
                    <div class="setting-group">
                        <label>Peak Labels</label>
                        <select id="peak-label-mode" onchange="renderSpectra()">
                            <option value="none">None</option>
                            <option value="top5">Top 5</option>
                            <option value="top10">Top 10</option>
                            <option value="all">All Major</option>
                        </select>
                    </div>
                    <div class="setting-group">
                        <label class="checkbox-inline">
                            <input type="checkbox" id="invert-x-axis" checked onchange="renderSpectra()">
                            Invert X-Axis (IR)
                        </label>
                    </div>
                </div>
                
                <div class="setting-group" style="margin-top: 10px; width: 100%;">
                    <label><strong>Select Molecules to Compare:</strong></label>
                    <div id="spectra-multiselect" class="multiselect-grid">
                        <!-- Checkboxes injected by JS -->
                    </div>
                </div>
            </details>
            
            <div class="tabs">
                <div class="tab active" onclick="showSpectraTab('ir')">🔴 IR Spectrum</div>
                <div class="tab" onclick="showSpectraTab('raman')">🟢 Raman Spectrum</div>
                <div class="tab" onclick="showSpectraTab('uvvis')">🟣 UV-Vis (TDDFT)</div>
                <div class="tab" onclick="showSpectraTab('correlation')">🔗 IR-Raman</div>
            </div>
            <div id="spectra-ir" class="tab-content active">
                <div id="ir-chart" class="chart-container"></div>
            </div>
            <div id="spectra-raman" class="tab-content">
                <div id="raman-chart" class="chart-container"></div>
            </div>
            <div id="spectra-uvvis" class="tab-content">
                <!-- UV-Vis specific settings -->
                <details class="settings-panel">
                    <summary>UV-Vis Settings</summary>
                    <div class="settings-row">
                        <div class="setting-group">
                            <label>Wavelength Range (nm)</label>
                            <div class="dual-range">
                                <input type="number" id="uvvis-wl-min" value="200" min="100" max="900" onchange="renderSpectra()">
                                <span>to</span>
                                <input type="number" id="uvvis-wl-max" value="700" min="100" max="1000" onchange="renderSpectra()">
                            </div>
                        </div>
                        <div class="setting-group">
                            <label>FWHM Broadening (nm)</label>
                            <input type="range" id="uvvis-fwhm" min="5" max="80" value="20" onchange="renderSpectra()">
                            <span class="range-value" id="uvvis-fwhm-val">20 nm</span>
                        </div>
                        <div class="setting-group">
                            <label class="checkbox-inline">
                                <input type="checkbox" id="uvvis-show-sticks" checked onchange="renderSpectra()">
                                Show Stick Spectrum
                            </label>
                        </div>
                    </div>
                </details>
                <div id="uvvis-chart" class="chart-container"></div>
            </div>
            <div id="spectra-correlation" class="tab-content">
                <!-- IR-Raman Correlation Settings -->
                <details class="settings-panel" open>
                    <summary>Correlation Settings</summary>
                    <div class="settings-row">
                        <div class="setting-group">
                            <label>Select Molecule</label>
                            <select id="corr-mol-select" onchange="renderCorrelation()">
                                <!-- Options injected by JS -->
                            </select>
                        </div>
                        <div class="setting-group">
                            <label>Max Pairing Distance (cm⁻¹)</label>
                            <input type="range" id="corr-max-delta" min="10" max="100" value="40" onchange="renderCorrelation()">
                            <span class="range-value" id="corr-delta-val">40</span>
                        </div>
                        <div class="setting-group">
                            <label>Peak Threshold (%)</label>
                            <input type="range" id="corr-threshold" min="1" max="50" value="5" onchange="renderCorrelation()">
                            <span class="range-value" id="corr-thresh-val">5%</span>
                        </div>
                    </div>
                    <div class="settings-row">
                        <div class="setting-group">
                            <label class="checkbox-inline">
                                <input type="checkbox" id="corr-connectors" checked onchange="renderCorrelation()">
                                Show Connector Lines
                            </label>
                        </div>
                        <div class="setting-group">
                            <label class="checkbox-inline">
                                <input type="checkbox" id="corr-invert-x" checked onchange="renderCorrelation()">
                                Invert X-Axis
                            </label>
                        </div>
                        <div class="setting-group">
                            <label>Peak Labels</label>
                            <select id="corr-peak-mode" onchange="renderCorrelation()">
                                <option value="paired">Paired Only</option>
                                <option value="all">All Peaks</option>
                                <option value="none">None</option>
                            </select>
                        </div>
                    </div>
                </details>
                <div id="correlation-chart" class="chart-container"></div>
                <div id="correlation-stats" style="display: flex; gap: 20px; justify-content: center; margin-top: 15px;">
                    <div class="stat-card" style="padding: 15px 25px;"><div class="stat-value" id="corr-ir-peaks">0</div><div class="stat-label">IR Peaks</div></div>
                    <div class="stat-card" style="padding: 15px 25px;"><div class="stat-value" id="corr-raman-peaks">0</div><div class="stat-label">Raman Peaks</div></div>
                    <div class="stat-card" style="padding: 15px 25px;"><div class="stat-value" id="corr-paired">0</div><div class="stat-label">Paired</div></div>
                </div>
            </div>
        </section>
'''

    
    if include_energy:
        html += '''
        <!-- 5. Energy Analysis -->
        <section id="energy-analysis" class="section">
            <h2>⚡ 5. Energy Analysis</h2>
            
            <p>
                Comparative energy analysis across molecular systems. View absolute or relative energies 
                and HOMO-LUMO gap diagrams.
            </p>
            
            <details class="settings-panel" open>
                <summary>Energy Comparison Settings</summary>
                <div class="settings-row">
                    <div class="setting-group">
                        <label class="checkbox-inline">
                            <input type="checkbox" id="energy-relative" checked onchange="renderEnergyChart()">
                            Show Relative Energy (kcal/mol)
                        </label>
                    </div>
                    <div class="setting-group">
                        <label>Energy Type</label>
                        <select id="energy-type" onchange="renderEnergyChart()">
                            <option value="auto">Auto (Gibbs or SP)</option>
                            <option value="gibbs">Gibbs Only</option>
                            <option value="sp">Single Point Only</option>
                        </select>
                    </div>
                </div>
                <div class="setting-group" style="width: 100%; margin-top: 10px;">
                    <label><strong>Select Molecules to Compare:</strong></label>
                    <div id="energy-multiselect" class="multiselect-grid">
                        <!-- Checkboxes injected by JS -->
                    </div>
                </div>
            </details>
            
            <div class="tabs">
                <div class="tab active" onclick="showEnergyTab('comparison')">📊 Energy Comparison</div>
                <div class="tab" onclick="showEnergyTab('homolumo')">🔋 HOMO-LUMO Diagram</div>
            </div>
            
            <div id="energy-comparison" class="tab-content active">
                <div id="energy-chart" class="chart-container"></div>
                <div id="energy-reference" style="text-align: center; margin-top: 10px; font-size: 13px; color: var(--text-secondary);"></div>
            </div>
            
            <div id="energy-homolumo" class="tab-content">
                <div id="homolumo-chart" class="chart-container"></div>
                <p style="font-size: 12px; color: var(--text-secondary); margin-top: 10px; text-align: center;">
                    Solid lines = HOMO (occupied), Dashed lines = LUMO (virtual). Gap values shown between levels.
                </p>
            </div>
        </section>
'''

    
    # Data table
    html += f'''        
        <!-- 7. Data Appendix -->
        <section id="data-appendix" class="section">
            <h2>📊 7. Complete Data Appendix</h2>
            
            <p>
                Complete tabulated data for all calculations including molecular identities, electronic states, 
                energies, and orbital energy levels.
            </p>
            
            <div style="overflow-x: auto;">
                <table class="data-table">
                    <thead>
                        <tr>
                            <th>Molecule</th>
                            <th>State</th>
                            <th>SMILES</th>
                            <th>Energy (Eh)</th>
                            <th>HOMO (eV)</th>
                            <th>LUMO (eV)</th>
                            <th>Gap (eV)</th>
                            <th>Method</th>
                        </tr>
                    </thead>
                    <tbody>
'''
    
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "N/A")
        state = row.get("optimized_state", "N/A")
        smiles = row.get("smiles", "N/A")
        energy = row.get("gibbs_Eh") or row.get("single_point_Eh")
        homo = row.get("homo_energy")
        lumo = row.get("lumo_energy")
        gap = row.get("homo_lumo_gap")
        method = row.get("method_id", "N/A")
        
        smiles_display = str(smiles)[:60] + "..." if smiles and len(str(smiles)) > 60 else smiles
        
        html += f'''                        <tr>
                            <td><strong>{mol_id}</strong></td>
                            <td>{state}</td>
                            <td class="smiles-cell" title="{smiles}">{smiles_display if smiles and str(smiles) != "nan" else "N/A"}</td>
                            <td>{f'{energy:.6f}' if energy else 'N/A'}</td>
                            <td>{f'{homo:.3f}' if homo else 'N/A'}</td>
                            <td>{f'{lumo:.3f}' if lumo else 'N/A'}</td>
                            <td>{f'{gap:.3f}' if gap else 'N/A'}</td>
                            <td>{method}</td>
                        </tr>
'''
    
    html += f'''                    </tbody>
                </table>
            </div>
        </section>
        
        <footer>
            <p>Generated by <strong>ORCA Visualization Platform</strong></p>
            <p>Report created: {pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
            <p style="margin-top: 10px; font-size: 12px;">
                This report was automatically generated from ORCA quantum chemistry calculation results.
            </p>
        </footer>
    </div>
    
    <script>
        // Embedded molecule data
        const molecules = {molecules_json};
        
        let viewer = null;
        let spinning = false;
        let spinInterval = null;
        let currentStyle = 'ballstick';
        let currentMolIndex = 0;
        
        // Initialize on load
        document.addEventListener('DOMContentLoaded', function() {{
            // Initialize viewer and molecule data
            if (molecules.length > 0) {{
                initViewer();
            }} else {{
                console.warn("No molecule data found!");
            }}
        }});
        
        function initViewer() {{
            console.log("InitViewer: Starting...");
            console.log("InitViewer: Molecule count =", molecules.length);
            
            if (molecules.length === 0) {{
                console.error("InitViewer: No molecules found!");
                return;
            }}
            
            console.log("InitViewer: First molecule =", molecules[0]);
            
            const viewerDiv = document.getElementById('viewer-3d');
            if (!viewerDiv) {{
                console.warn("InitViewer: viewer-3d element not found");
            }} else {{
                viewer = $3Dmol.createViewer('viewer-3d', {{
                    backgroundColor: '{"#1e1e2e" if dark_theme else "#f0f0f5"}'
                }});
                console.log("InitViewer: 3Dmol viewer created");
            }}
            
            // Generate checkboxes
            generateCheckboxes('spectra-multiselect', renderSpectra);
            generateCheckboxes('energy-multiselect', renderEnergyChart);
            console.log("InitViewer: Checkboxes generated");
            
            // Populate molecule dropdowns
            const orbSelect = document.getElementById('orbital-mol-select');
            const corrSelect = document.getElementById('corr-mol-select');
            molecules.forEach((mol, idx) => {{
                if (orbSelect) {{
                    const opt = document.createElement('option');
                    opt.value = idx;
                    opt.textContent = mol.label;
                    orbSelect.appendChild(opt);
                }}
                if (corrSelect) {{
                    const opt = document.createElement('option');
                    opt.value = idx;
                    opt.textContent = mol.label;
                    corrSelect.appendChild(opt);
                }}
            }});
            console.log("InitViewer: Molecule dropdowns populated");
            
            loadMolecule(0);
            console.log("InitViewer: Complete");
        }}
        
        function loadMolecule(index) {{
            console.log("loadMolecule: Starting for index", index);
            
            if (!molecules[index]) {{
                console.error("loadMolecule: No molecule at index", index);
                return;
            }}
            
            currentMolIndex = index;
            const mol = molecules[index];
            console.log("loadMolecule: Processing", mol.label);
            
            // 3D Viewer (may fail if $3Dmol not loaded)
            try {{
                if (viewer && mol.xyz) {{
                    viewer.removeAllModels();
                    viewer.addModel(mol.xyz, 'xyz');
                    applyStyle(currentStyle);
                    viewer.zoomTo();
                    viewer.render();
                    console.log("loadMolecule: 3D viewer updated");
                }} else if (!viewer) {{
                    console.warn("loadMolecule: No 3D viewer available");
                }} else {{
                    console.warn("loadMolecule: No XYZ data for", mol.label);
                }}
            }} catch (e) {{
                console.error("loadMolecule: 3D viewer error:", e);
            }}
            
            // 2D Image
            try {{
                if (mol.svg) {{
                    const img = document.getElementById('img-2d');
                    if (img) {{
                        img.src = 'data:image/svg+xml;base64,' + mol.svg;
                        img.style.display = 'block';
                    }}
                    const fallback = document.getElementById('no-2d-msg');
                    if (fallback) fallback.style.display = 'none';
                }} else {{
                    const img = document.getElementById('img-2d');
                    if (img) img.style.display = 'none';
                    const fallback = document.getElementById('no-2d-msg');
                    if (fallback) fallback.style.display = 'block';
                }}
            }} catch (e) {{
                console.error("loadMolecule: 2D image error:", e);
            }}
            
            // Update molecule info panel
            try {{
                updateMolInfo(index);
                console.log("loadMolecule: Mol info updated");
            }} catch (e) {{
                console.error("loadMolecule: updateMolInfo error:", e);
            }}
            
            // Ensure current molecule is selected in comparison views
            ['spectra-multiselect', 'energy-multiselect'].forEach(id => {{
                 const div = document.getElementById(id);
                 if(div) {{
                     const box = div.querySelector(`input[value="${{index}}"]`);
                     if(box) box.checked = true;
                 }}
            }});
            
            // Render all charts (independent of 3D viewer)
            try {{
                renderSpectra();
                console.log("loadMolecule: Spectra rendered");
            }} catch (e) {{
                console.error("loadMolecule: renderSpectra error:", e);
            }}
            
            try {{
                renderOrbitals(index);
                console.log("loadMolecule: Orbitals rendered");
            }} catch (e) {{
                console.error("loadMolecule: renderOrbitals error:", e);
            }}
            
            try {{
                renderEnergyChart();
                console.log("loadMolecule: Energy chart rendered");
            }} catch (e) {{
                console.error("loadMolecule: renderEnergyChart error:", e);
            }}
            
            console.log("loadMolecule: Complete");
        }}
        
        function switchMolecule() {{
            const select = document.getElementById('molecule-select');
            loadMolecule(parseInt(select.value));
        }}
        
        function setStyle(style) {{
            currentStyle = style;
            applyStyle(style);
            
            // Update button states
            document.querySelectorAll('.controls button').forEach(b => b.classList.remove('active'));
            event.target.classList.add('active');
        }}
        
        function applyStyle(style) {{
            if (!viewer) return;
            
            let spec = {{}};
            let common = {{}};
            
            // Apply color scheme if not default
            if (colorScheme !== 'default') {{
                common.colorscheme = colorScheme;
            }}
            
            if (style === 'stick') {{
                spec = {{stick: {{radius: bondRadius * 1.5, ...common}}}};
            }} else if (style === 'sphere') {{
                spec = {{sphere: {{scale: atomScale * 3.0, ...common}}}};
            }} else {{
                // Ball and Stick
                spec = {{
                    stick: {{radius: bondRadius, ...common}},
                    sphere: {{scale: atomScale, ...common}}
                }};
            }}
            
            viewer.setStyle({{}}, spec);
            viewer.render();
        }}
        
        function toggleSpin() {{
            spinning = !spinning;
            const btn = document.getElementById('btn-spin');
            if (spinning) {{
                spinInterval = setInterval(() => {{
                    viewer.rotate(1, {{x: 0, y: 1, z: 0}});
                    viewer.render();
                }}, 40);
                if (btn) btn.classList.add('active');
            }} else {{
                clearInterval(spinInterval);
                if (btn) btn.classList.remove('active');
            }}
        }}
        
        // Settings handlers
        let showLabels = false;
        let bondRadius = 0.12;
        let atomScale = 0.25;
        let colorScheme = 'Jmol';
        
        function toggleLabels() {{
            showLabels = !showLabels;
            const btn = document.getElementById('btn-labels');
            if (showLabels) {{
                // Add labels to non-hydrogen atoms
                const mol = molecules[currentMolIndex];
                if (mol && mol.coords) {{
                    mol.coords.forEach((atom, i) => {{
                        if (atom.el !== 'H') {{
                            viewer.addLabel(atom.el + (i+1), {{
                                position: {{x: atom.x, y: atom.y, z: atom.z}},
                                fontSize: 12,
                                fontColor: 'white',
                                backgroundOpacity: 0.6,
                                backgroundColor: '#333'
                            }});
                        }}
                    }});
                }}
                if (btn) btn.classList.add('active');
            }} else {{
                viewer.removeAllLabels();
                if (btn) btn.classList.remove('active');
            }}
            viewer.render();
        }}
        
        function updateBondRadius(val) {{
            bondRadius = parseFloat(val);
            document.getElementById('bond-radius-val').textContent = val;
            applyStyle(currentStyle);
        }}
        
        function updateAtomScale(val) {{
            atomScale = parseFloat(val);
            document.getElementById('atom-scale-val').textContent = val;
            applyStyle(currentStyle);
        }}
        
        function updateColorScheme(scheme) {{
            colorScheme = scheme;
            applyStyle(currentStyle);
        }}
        
        function updateBackground(color) {{
            if (viewer) {{
                viewer.setBackgroundColor(color);
                viewer.render();
            }}
        }}

        function showMolTab(mode) {{
            // Hide all mol view contents
            ['mol-view-3d', 'mol-view-2d', 'mol-view-props'].forEach(id => {{
                const el = document.getElementById(id);
                if (el) el.classList.remove('active');
            }});
            
            // Deactivate all mol view tabs
            document.querySelectorAll('.viewer-container .tabs .tab').forEach(t => t.classList.remove('active'));
            
            // Activate clicked tab
            if (event && event.target) event.target.classList.add('active');
            
            // Show content
            const targetEl = document.getElementById('mol-view-' + mode);
            if (targetEl) targetEl.classList.add('active');
        }}
        
        function showSpectraTab(type) {{
            // Scope to spectra section only
            const section = document.getElementById('spectra-analysis');
            if (!section) return;
            
            section.querySelectorAll('.tabs .tab').forEach(t => t.classList.remove('active'));
            ['spectra-ir', 'spectra-raman', 'spectra-uvvis', 'spectra-correlation'].forEach(id => {{
                const el = document.getElementById(id);
                if (el) el.classList.remove('active');
            }});
            
            if (event && event.target) event.target.classList.add('active');
            const targetEl = document.getElementById('spectra-' + type);
            if (targetEl) targetEl.classList.add('active');
            
            // If switching to correlation, trigger render
            if (type === 'correlation') {{
                renderCorrelation();
            }}
        }}
        
        // Orbital variables
        let orbitalRange = 6;
        
        function updateOrbitalRange(val) {{
            orbitalRange = parseInt(val);
            document.getElementById('orb-range-val').textContent = val;
            renderOrbitals(currentMolIndex);
        }}
        
        function renderOrbitals(index) {{
            const div = document.getElementById('orbital-chart');
            if (!div) return;
            
            const mol = molecules[index];
            if (!mol || !mol.orbitals || mol.orbitals.length === 0) {{
                div.innerHTML = '<div style="text-align:center; padding: 20px;">No orbital data available</div>';
                return;
            }}
            
            // Sort by energy
            let orbs = [...mol.orbitals].sort((a, b) => a.energy - b.energy);
            const homoval = mol.homo || (orbs.length > 0 ? orbs[orbs.length/2].energy : 0);
            
            // Filter to N closest to HOMO
            // Calculate distance to HOMO
            orbs.forEach(o => o.dist = Math.abs(o.energy - homoval));
            orbs.sort((a, b) => a.dist - b.dist);
            
            // Take top N*2
            let subset = orbs.slice(0, orbitalRange * 2);
            // Sort back by energy for drawing
            subset.sort((a, b) => a.energy - b.energy);
            
            const shapes = [];
            const annotations = [];
            
            subset.forEach(o => {{
                const isOcc = (o.occ !== undefined && o.occ > 0.1) || (mol.homo !== undefined && o.energy <= mol.homo + 0.001);
                const color = isOcc ? '#3366cc' : '#cccccc';
                const labelColor = isOcc ? '#3366cc' : '#999999';
                
                // Draw line
                shapes.push({{
                    type: 'line',
                    x0: 0.2, x1: 0.8,
                    y0: o.energy, y1: o.energy,
                    line: {{color: color, width: 3}}
                }});
                
                // Label (Index/Occ)
                let text = "";
                if (o.occ !== undefined) text += `Occ: ${{o.occ.toFixed(1)}}`;
                
                annotations.push({{
                    x: 0.82, y: o.energy,
                    text: `${{o.energy.toFixed(2)}} eV`,
                    showarrow: false,
                    xanchor: 'left',
                    font: {{size: 10, color: labelColor}}
                }});
            }});
            
            const layout = {{
                title: `Orbital Energy Levels (${{mol.label}})`,
                xaxis: {{showgrid: false, zeroline: false, showticklabels: false, range: [0, 1]}},
                yaxis: {{title: 'Energy (eV)'}},
                shapes: shapes,
                annotations: annotations,
                paper_bgcolor: '{"#2a2a3e" if dark_theme else "#fff"}',
                plot_bgcolor: '{"#2a2a3e" if dark_theme else "#fff"}',
                font: {{color: '{text_primary}'}},
                margin: {{l: 60, r: 100, t: 50, b: 30}}
            }};
            
            Plotly.newPlot('orbital-chart', [], layout, {{responsive: true}});
        }}
        
        // Spectra variables
        let spectraFWHM = 10;
        let spectraMode = 'line';
        
        function updateSpectraBroadening(val) {{
            spectraFWHM = parseFloat(val);
            document.getElementById('fwhm-val').textContent = val + " cm⁻¹";
            renderSpectra(currentMolIndex);
        }}
        
        function updateSpectraMode(val) {{
            spectraMode = val;
            renderSpectra();
        }}
        
        function gaussianBroadening(peaks, fwhm, xMin, xMax, nPoints=2000) {{
            const x = [];
            const y = [];
            const step = (xMax - xMin) / (nPoints - 1);
            const alpha = 4 * Math.log(2) / (fwhm * fwhm);
            
            for (let i = 0; i < nPoints; i++) {{
                const val = xMin + i * step;
                x.push(val);
                let intensity = 0;
                
                for (const p of peaks) {{
                    const dist = (val - p.freq);
                    if (Math.abs(dist) > 3 * fwhm) continue;
                    intensity += p.intensity * Math.exp(-alpha * dist * dist);
                }}
                y.push(intensity);
            }}
            return {{x, y}};
        }}
        
        // Multi-select helpers
        function generateCheckboxes(containerId, onChangeHandler) {{
            const div = document.getElementById(containerId);
            if (!div) return;
            
            div.innerHTML = "";
            molecules.forEach((mol, idx) => {{
                const label = document.createElement('label');
                label.className = 'multiselect-item';
                
                const input = document.createElement('input');
                input.type = 'checkbox';
                input.value = idx;
                input.checked = (idx === currentMolIndex); // Default check current
                input.onchange = () => onChangeHandler();
                
                label.appendChild(input);
                label.appendChild(document.createTextNode(mol.label));
                div.appendChild(label);
            }});
        }}
        
        function getSelectedIndices(containerId) {{
            const div = document.getElementById(containerId);
            if (!div) return [currentMolIndex];
            
            const checkboxes = div.querySelectorAll('input[type="checkbox"]:checked');
            if (checkboxes.length === 0) return [currentMolIndex]; // Fallback to current
            
            return Array.from(checkboxes).map(cb => parseInt(cb.value));
        }}
        
        function renderSpectra() {{
            const indices = getSelectedIndices('spectra-multiselect');
            const colors = ['#00ff88', '#ff6b6b', '#6b88ff', '#ffd93d', '#ff9f43', '#a55eea', '#2bcbba', '#fd9644'];
            
            // IR Chart
            const irDiv = document.getElementById('ir-chart');
            if (irDiv) {{
                const traces = [];
                indices.forEach((idx, i) => {{
                    const mol = molecules[idx];
                    if (!mol.ir || mol.ir.length === 0) return;
                    
                    const color = colors[i % colors.length];
                    let trace = {{}};
                    
                    if (spectraMode === 'line') {{
                       const freqs = mol.ir.map(p => p.freq);
                       const minFreq = Math.min(...freqs) - 200;
                       const maxFreq = Math.max(...freqs) + 200;
                       const broadened = gaussianBroadening(mol.ir, spectraFWHM, Math.max(0, minFreq), maxFreq);
                       
                       trace = {{
                           x: broadened.x, y: broadened.y,
                           type: 'scatter', mode: 'lines',
                           line: {{color: color, width: 2}},
                           name: mol.label
                       }};
                       // Only fill if single trace
                       if (indices.length === 1) {{
                           trace.fill = 'tozeroy';
                           trace.fillcolor = color + '1A'; // 10% opacity
                       }}
                    }} else {{
                        const x = [];
                        const y = [];
                        mol.ir.forEach(p => {{
                            x.push(p.freq, p.freq, null);
                            y.push(0, p.intensity, null);
                        }});
                        trace = {{
                            x: x, y: y,
                            type: 'scatter', mode: 'lines',
                            line: {{color: color, width: 2}},
                            name: mol.label
                        }};
                    }}
                    traces.push(trace);
                }});
                
                Plotly.newPlot('ir-chart', traces, {{
                    title: 'IR Spectrum',
                    xaxis: {{title: 'Wavenumber (cm⁻¹)', autorange: document.getElementById('invert-x-axis')?.checked ? 'reversed' : true}},
                    yaxis: {{title: 'Intensity'}},
                    paper_bgcolor: '{"#2a2a3e" if dark_theme else "#fff"}',
                    plot_bgcolor: '{"#2a2a3e" if dark_theme else "#fff"}',
                    font: {{color: '{text_primary}'}},
                    showlegend: indices.length > 1,
                    margin: {{t: 30, b: 30, l: 50, r: 20}}
                }}, {{responsive: true}});
            }}
            
            // Raman Chart
            const ramanDiv = document.getElementById('raman-chart');
            if (ramanDiv) {{
                const traces = [];
                indices.forEach((idx, i) => {{
                    const mol = molecules[idx];
                    if (!mol.raman || mol.raman.length === 0) return;
                    
                    const color = colors[i % colors.length];
                    let trace = {{}};
                    
                    if (spectraMode === 'line') {{
                        const freqs = mol.raman.map(p => p.freq);
                        const minFreq = Math.min(...freqs) - 200;
                        const maxFreq = Math.max(...freqs) + 200;
                        const peaks = mol.raman.map(p => ({{freq: p.freq, intensity: p.activity}}));
                        const broadened = gaussianBroadening(peaks, spectraFWHM, Math.max(0, minFreq), maxFreq);
                        
                        trace = {{
                            x: broadened.x, y: broadened.y,
                            type: 'scatter', mode: 'lines',
                            line: {{color: color, width: 2}},
                            name: mol.label
                        }};
                        if (indices.length === 1) {{
                            trace.fill = 'tozeroy';
                            trace.fillcolor = color + '1A';
                        }}
                    }} else {{
                        const x = [];
                        const y = [];
                        mol.raman.forEach(p => {{
                           x.push(p.freq, p.freq, null);
                           y.push(0, p.activity, null);
                        }});
                        trace = {{
                           x: x, y: y,
                           type: 'scatter', mode: 'lines',
                           line: {{color: color, width: 2}},
                           name: mol.label
                        }};
                    }}
                    traces.push(trace);
                }});
                
                Plotly.newPlot('raman-chart', traces, {{
                    title: 'Raman Spectrum',
                    xaxis: {{title: 'Shift (cm⁻¹)'}},
                    yaxis: {{title: 'Intensity'}},
                    paper_bgcolor: '{"#2a2a3e" if dark_theme else "#fff"}',
                    plot_bgcolor: '{"#2a2a3e" if dark_theme else "#fff"}',
                    font: {{color: '{text_primary}'}},
                    showlegend: indices.length > 1,
                    margin: {{t: 30, b: 30, l: 50, r: 20}}
                }}, {{responsive: true}});
            }}
            
             // UV-Vis Chart
            const uvDiv = document.getElementById('uvvis-chart');
            if (uvDiv) {{
                 const traces = [];
                 const uvFwhm = parseFloat(document.getElementById('uvvis-fwhm')?.value) || 20;
                 const wlMin = parseFloat(document.getElementById('uvvis-wl-min')?.value) || 200;
                 const wlMax = parseFloat(document.getElementById('uvvis-wl-max')?.value) || 700;
                 const showSticks = document.getElementById('uvvis-show-sticks')?.checked ?? true;
                 
                 // Update FWHM display
                 const fwhmVal = document.getElementById('uvvis-fwhm-val');
                 if (fwhmVal) fwhmVal.textContent = uvFwhm + ' nm';
                 
                 indices.forEach((idx, i) => {{
                    const mol = molecules[idx];
                    if (!mol.tddft || mol.tddft.length === 0) return;
                    
                    const color = colors[i % colors.length];
                    
                    // Convert TDDFT states to peaks for broadening
                    const peaks = mol.tddft.filter(s => s.nm && s.f).map(s => ({{freq: s.nm, intensity: s.f}}));
                    
                    if (peaks.length === 0) return;
                    
                    // Gaussian broadening on wavelength
                    const broadened = gaussianBroadening(peaks, uvFwhm, wlMin, wlMax);
                    
                    traces.push({{
                        x: broadened.x, y: broadened.y,
                        type: 'scatter', mode: 'lines',
                        line: {{color: color, width: 2}},
                        fill: indices.length === 1 ? 'tozeroy' : 'none',
                        fillcolor: color + '1A',
                        name: mol.label
                    }});
                    
                    // Add stick spectrum if enabled
                    if (showSticks) {{
                        const stickX = [];
                        const stickY = [];
                        peaks.forEach(p => {{
                            stickX.push(p.freq, p.freq, null);
                            stickY.push(0, p.intensity, null);
                        }});
                        traces.push({{
                            x: stickX, y: stickY,
                            type: 'scatter', mode: 'lines',
                            line: {{color: color, width: 1.5}},
                            showlegend: false,
                            hoverinfo: 'skip'
                        }});
                    }}
                 }});

                 Plotly.newPlot('uvvis-chart', traces, {{
                    title: 'UV-Vis Absorption Spectrum',
                    xaxis: {{title: 'Wavelength (nm)', range: [wlMin, wlMax]}},
                    yaxis: {{title: 'Oscillator Strength'}},
                    paper_bgcolor: '{"#2a2a3e" if dark_theme else "#fff"}',
                    plot_bgcolor: '{"#2a2a3e" if dark_theme else "#fff"}',
                    font: {{color: '{text_primary}'}},
                    showlegend: indices.length > 1,
                    margin: {{t: 30, b: 30, l: 50, r: 20}}
                }}, {{responsive: true}});
            }}
        }}
        
        

        
        function renderEnergyChart() {{
            const indices = getSelectedIndices('energy-multiselect');
            const showRelative = document.getElementById('energy-relative')?.checked ?? false;
            const energyType = document.getElementById('energy-type')?.value || 'auto';
            
            // Get energy based on type
            function getEnergy(mol) {{
                if (energyType === 'gibbs') return mol.gibbs;
                if (energyType === 'sp') return mol.single_point;
                // auto: prefer gibbs, fallback to sp
                return mol.gibbs || mol.single_point || mol.energy;
            }}
            
            const x = [];
            let y = [];
            const markerColors = [];
            const text = [];
            
            indices.forEach(idx => {{
                const mol = molecules[idx];
                x.push(mol.label);
                y.push(getEnergy(mol) || 0);
                markerColors.push(idx === currentMolIndex ? '#667eea' : '#a3bffa');
            }});
            
            // Apply relative energy calculation
            let yLabel = 'Energy (Eh)';
            let refLabel = '';
            
            if (showRelative && y.length > 0 && y.some(v => v !== 0)) {{
                const minE = Math.min(...y.filter(v => v !== 0));
                const refIdx = y.indexOf(minE);
                refLabel = x[refIdx];
                
                // Convert to relative kcal/mol
                y = y.map(e => {{
                    if (e === 0) return 0;
                    return (e - minE) * 627.509; // Eh to kcal/mol
                }});
                yLabel = 'Relative Energy (kcal/mol)';
            }}
            
            // Generate text labels
            y.forEach((energy, i) => {{
                if (showRelative) {{
                    text.push(`${{energy.toFixed(2)}} kcal/mol`);
                }} else {{
                    text.push(`${{energy.toFixed(5)}} Eh`);
                }}
            }});
            
            const layout = {{
                title: showRelative ? 'Relative Energy Comparison' : 'Energy Comparison',
                xaxis: {{title: 'Molecule'}},
                yaxis: {{title: yLabel}},
                paper_bgcolor: '{"#2a2a3e" if dark_theme else "#fff"}',
                plot_bgcolor: '{"#2a2a3e" if dark_theme else "#fff"}',
                font: {{color: '{text_primary}'}},
                margin: {{l: 80, r: 30, t: 50, b: 80}}
            }};
            
            Plotly.newPlot('energy-chart', [{{
                x: x, y: y,
                type: 'bar',
                marker: {{color: markerColors}},
                text: text,
                textposition: 'outside',
                hoverinfo: 'x+text'
            }}], layout, {{responsive: true}});
            
            // Update reference label
            const refDiv = document.getElementById('energy-reference');
            if (refDiv) {{
                if (showRelative && refLabel) {{
                    refDiv.innerHTML = `<strong>Reference:</strong> ${{refLabel}} (set to 0.00 kcal/mol)`;
                }} else {{
                    refDiv.innerHTML = '';
                }}
            }}
        }}
        
        function updateMolInfo(index) {{
            const mol = molecules[index];
            if (!mol) return;
            
            let html = `
                <div class="stat-card"><div class="stat-value">${{mol.atoms || 0}}</div><div class="stat-label">Atoms</div></div>
                <div class="stat-card"><div class="stat-value">${{mol.energy ? mol.energy.toFixed(2) : 'N/A'}}</div><div class="stat-label">Energy (Eh)</div></div>
                <div class="stat-card"><div class="stat-value">${{mol.homo ? mol.homo.toFixed(2) : 'N/A'}}</div><div class="stat-label">HOMO (eV)</div></div>
                <div class="stat-card"><div class="stat-value">${{mol.lumo ? mol.lumo.toFixed(2) : 'N/A'}}</div><div class="stat-label">LUMO (eV)</div></div>
            `;
            
            // Add SMILES if available
            if (mol.smiles) {{
                html += `<div class="stat-card" style="grid-column: span 4; text-align: left;">
                    <strong>SMILES:</strong> <code style="font-size: 0.9em;">${{mol.smiles}}</code>
                </div>`;
            }}
            
            // Add method info
            if (mol.method || mol.functional) {{
                html += `<div class="stat-card" style="grid-column: span 4; text-align: left;">
                    <strong>Method:</strong> ${{mol.method || mol.functional || 'N/A'}} | <strong>Basis:</strong> ${{mol.basis_set || 'N/A'}}
                </div>`;
            }}
            
            document.getElementById('mol-info').innerHTML = html;
            
            // Populate property tables
            populatePropertyTables(mol);
        }}
        
        function populatePropertyTables(mol) {{
            // Coordinates table
            const coordsTbody = document.getElementById('coords-tbody');
            if (coordsTbody && mol.coords) {{
                let html = '';
                mol.coords.forEach((atom, i) => {{
                    html += `<tr><td>${{i+1}}</td><td>${{atom.el}}</td><td>${{atom.x.toFixed(4)}}</td><td>${{atom.y.toFixed(4)}}</td><td>${{atom.z.toFixed(4)}}</td></tr>`;
                }});
                coordsTbody.innerHTML = html || '<tr><td colspan="5">No coordinate data</td></tr>';
            }}
            
            // Properties table
            const propsTbody = document.getElementById('props-tbody');
            if (propsTbody) {{
                let html = '';
                const props = [
                    ['Molecule ID', mol.label],
                    ['State', mol.state || 'N/A'],
                    ['SMILES', mol.smiles || 'N/A'],
                    ['Energy (Eh)', mol.energy ? mol.energy.toFixed(6) : 'N/A'],
                    ['HOMO (eV)', mol.homo ? mol.homo.toFixed(4) : 'N/A'],
                    ['LUMO (eV)', mol.lumo ? mol.lumo.toFixed(4) : 'N/A'],
                    ['Gap (eV)', mol.gap ? mol.gap.toFixed(4) : 'N/A'],
                    ['Charge', mol.charge !== undefined ? mol.charge : 'N/A'],
                    ['Multiplicity', mol.multiplicity || 'N/A']
                ];
                props.forEach(([k, v]) => {{
                    html += `<tr><td><strong>${{k}}</strong></td><td>${{v}}</td></tr>`;
                }});
                propsTbody.innerHTML = html;
            }}
            
            // Mulliken charges table
            const mullikenTbody = document.getElementById('mulliken-tbody');
            if (mullikenTbody && mol.mulliken) {{
                let html = '';
                mol.mulliken.forEach((item, i) => {{
                    const charge = item.charge !== undefined ? item.charge.toFixed(4) : 'N/A';
                    const spin = item.spin !== undefined ? item.spin.toFixed(4) : 'N/A';
                    html += `<tr><td>${{i+1}}</td><td>${{item.atom || item.el || 'X'}}</td><td>${{charge}}</td><td>${{spin}}</td></tr>`;
                }});
                mullikenTbody.innerHTML = html || '<tr><td colspan="4">No Mulliken data</td></tr>';
            }} else if (mullikenTbody) {{
                mullikenTbody.innerHTML = '<tr><td colspan="4">No Mulliken data</td></tr>';
            }}
            
            // Method table
            const methodTbody = document.getElementById('method-tbody');
            if (methodTbody) {{
                let html = '';
                const methodProps = [
                    ['Functional', mol.functional || 'N/A'],
                    ['Basis Set', mol.basis_set || 'N/A'],
                    ['Dispersion', mol.dispersion || 'N/A'],
                    ['Solvent', mol.solvent || 'N/A'],
                    ['Method ID', mol.method_id || mol.method || 'N/A']
                ];
                methodProps.forEach(([k, v]) => {{
                    html += `<tr><td><strong>${{k}}</strong></td><td>${{v}}</td></tr>`;
                }});
                methodTbody.innerHTML = html;
            }}
        }}
        
        // Property sub-tab navigation
        function showPropSubTab(type) {{
            ['coords', 'properties', 'mulliken', 'method'].forEach(t => {{
                const el = document.getElementById('prop-' + t);
                if (el) el.classList.remove('active');
            }});
            document.querySelectorAll('.prop-tab').forEach(tab => tab.classList.remove('active'));
            
            const targetContent = document.getElementById('prop-' + type);
            if (targetContent) targetContent.classList.add('active');
            
            if (event && event.target) event.target.classList.add('active');
        }}
        
        // 2D Settings update
        function update2DSettings() {{
            const size = document.getElementById('mol2d-size')?.value || 450;
            const bg = document.getElementById('mol2d-bg')?.value || '#ffffff';
            const fontsize = document.getElementById('mol2d-fontsize')?.value || 14;
            
            // Update display values
            const sizeVal = document.getElementById('mol2d-size-val');
            if (sizeVal) sizeVal.textContent = size + 'px';
            const bonddVal = document.getElementById('mol2d-bondwidth-val');
            if (bonddVal) bonddVal.textContent = document.getElementById('mol2d-bondwidth')?.value || '2.0';
            const paddingVal = document.getElementById('mol2d-padding-val');
            if (paddingVal) paddingVal.textContent = document.getElementById('mol2d-padding')?.value || '0.05';
            const fontsizeVal = document.getElementById('mol2d-fontsize-val');
            if (fontsizeVal) fontsizeVal.textContent = fontsize;
            
            // Update container
            const container = document.getElementById('mol2d-container');
            if (container) {{
                container.style.background = bg;
            }}
        }}
        
        // Energy tab navigation
        function showEnergyTab(type) {{
            ['comparison', 'homolumo'].forEach(t => {{
                const el = document.getElementById('energy-' + t);
                if (el) el.classList.remove('active');
            }});
            
            const section = document.getElementById('energy-analysis');
            if (section) {{
                section.querySelectorAll('.tab').forEach(tab => tab.classList.remove('active'));
            }}
            
            const targetContent = document.getElementById('energy-' + type);
            if (targetContent) targetContent.classList.add('active');
            
            if (event && event.target) event.target.classList.add('active');
            
            // Render chart if switching to homolumo
            if (type === 'homolumo') {{
                renderHomoLumoChart();
            }}
        }}
        
        // HOMO-LUMO Diagram
        function renderHomoLumoChart() {{
            const div = document.getElementById('homolumo-chart');
            if (!div) return;
            
            const indices = getSelectedIndices('energy-multiselect');
            const colors = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880'];
            const traces = [];
            const annotations = [];
            
            indices.forEach((idx, i) => {{
                const mol = molecules[idx];
                if (!mol.homo && !mol.lumo) return;
                
                const color = colors[i % colors.length];
                const x = i;
                
                // HOMO line (solid)
                if (mol.homo !== undefined && mol.homo !== null) {{
                    traces.push({{
                        x: [x - 0.3, x + 0.3],
                        y: [mol.homo, mol.homo],
                        mode: 'lines',
                        line: {{color: color, width: 5}},
                        showlegend: false,
                        hovertemplate: `${{mol.label}}<br>HOMO: ${{mol.homo.toFixed(3)}} eV<extra></extra>`
                    }});
                }}
                
                // LUMO line (dashed)
                if (mol.lumo !== undefined && mol.lumo !== null) {{
                    traces.push({{
                        x: [x - 0.3, x + 0.3],
                        y: [mol.lumo, mol.lumo],
                        mode: 'lines',
                        line: {{color: color, width: 5, dash: 'dash'}},
                        showlegend: false,
                        hovertemplate: `${{mol.label}}<br>LUMO: ${{mol.lumo.toFixed(3)}} eV<extra></extra>`
                    }});
                }}
                
                // Connecting line with gap annotation
                if (mol.homo !== undefined && mol.lumo !== undefined) {{
                    traces.push({{
                        x: [x, x],
                        y: [mol.homo, mol.lumo],
                        mode: 'lines',
                        line: {{color: color, width: 1, dash: 'dot'}},
                        showlegend: false,
                        hoverinfo: 'skip'
                    }});
                    
                    const gap = mol.lumo - mol.homo;
                    annotations.push({{
                        x: x,
                        y: (mol.homo + mol.lumo) / 2,
                        text: `${{gap.toFixed(2)}} eV`,
                        showarrow: false,
                        font: {{size: 10}},
                        bgcolor: 'white',
                        borderpad: 2
                    }});
                }}
            }});
            
            const labels = indices.map(i => molecules[i].label);
            
            Plotly.newPlot(div, traces, {{
                title: 'HOMO-LUMO Energy Levels',
                xaxis: {{
                    tickvals: indices.map((_, i) => i),
                    ticktext: labels,
                    title: 'Molecule'
                }},
                yaxis: {{title: 'Energy (eV)'}},
                annotations: annotations,
                paper_bgcolor: '{"#2a2a3e" if dark_theme else "#fff"}',
                plot_bgcolor: '{"#2a2a3e" if dark_theme else "#fff"}',
                font: {{color: '{text_primary}'}},
                margin: {{l: 60, r: 30, t: 50, b: 80}}
            }}, {{responsive: true}});
        }}
        
        // IR-Raman Correlation
        function renderCorrelation() {{
            const div = document.getElementById('correlation-chart');
            if (!div) return;
            
            const molSelect = document.getElementById('corr-mol-select');
            if (!molSelect) return;
            
            // Populate molecule selector if empty
            if (molSelect.options.length === 0) {{
                molecules.forEach((mol, idx) => {{
                    const opt = document.createElement('option');
                    opt.value = idx;
                    opt.textContent = mol.label;
                    molSelect.appendChild(opt);
                }});
            }}
            
            const idx = parseInt(molSelect.value);
            const mol = molecules[idx];
            if (!mol) return;
            
            // Get settings
            const freqMin = parseInt(document.getElementById('freq-min')?.value) || 400;
            const freqMax = parseInt(document.getElementById('freq-max')?.value) || 4000;
            const maxDelta = parseInt(document.getElementById('corr-max-delta')?.value) || 40;
            const threshold = (parseInt(document.getElementById('corr-threshold')?.value) || 5) / 100;
            const showConnectors = document.getElementById('corr-connectors')?.checked ?? true;
            const invertX = document.getElementById('corr-invert-x')?.checked ?? true;
            
            // Update display values
            const deltaVal = document.getElementById('corr-delta-val');
            if (deltaVal) deltaVal.textContent = maxDelta;
            const threshVal = document.getElementById('corr-thresh-val');
            if (threshVal) threshVal.textContent = (threshold * 100) + '%';
            
            if (!mol.ir || mol.ir.length === 0 || !mol.raman || mol.raman.length === 0) {{
                div.innerHTML = '<div style="text-align:center; padding: 40px; color: var(--text-secondary);">No IR or Raman data available for this molecule</div>';
                return;
            }}
            
            // Apply broadening
            const fwhm = spectraFWHM;
            const irBroad = gaussianBroadening(mol.ir, fwhm, freqMin, freqMax);
            const ramanPeaks = mol.raman.map(p => ({{freq: p.freq, intensity: p.activity || p.intensity}}));
            const ramanBroad = gaussianBroadening(ramanPeaks, fwhm, freqMin, freqMax);
            
            // Normalize
            const irMax = Math.max(...irBroad.y);
            const ramanMax = Math.max(...ramanBroad.y);
            const irNorm = irBroad.y.map(v => v / irMax);
            const ramanNorm = ramanBroad.y.map(v => v / ramanMax);
            
            // Find peaks (simple local maxima)
            function findPeaks(y, thresh) {{
                const peaks = [];
                for (let i = 1; i < y.length - 1; i++) {{
                    if (y[i] > y[i-1] && y[i] > y[i+1] && y[i] > thresh) {{
                        peaks.push({{idx: i, val: y[i]}});
                    }}
                }}
                return peaks;
            }}
            
            const irPeaks = findPeaks(irNorm, threshold);
            const ramanPeaks2 = findPeaks(ramanNorm, threshold);
            
            // Auto-pair by frequency
            const pairs = [];
            const usedIr = new Set();
            const usedRaman = new Set();
            
            const candidates = [];
            irPeaks.forEach((ip, iIdx) => {{
                ramanPeaks2.forEach((rp, rIdx) => {{
                    const irFreq = irBroad.x[ip.idx];
                    const ramanFreq = ramanBroad.x[rp.idx];
                    const delta = Math.abs(irFreq - ramanFreq);
                    if (delta <= maxDelta) {{
                        candidates.push({{delta, iIdx, rIdx, irFreq, ramanFreq, irInt: ip.val, ramanInt: rp.val}});
                    }}
                }});
            }});
            
            candidates.sort((a, b) => a.delta - b.delta);
            candidates.forEach(c => {{
                if (!usedIr.has(c.iIdx) && !usedRaman.has(c.rIdx)) {{
                    pairs.push(c);
                    usedIr.add(c.iIdx);
                    usedRaman.add(c.rIdx);
                }}
            }});
            
            // Build traces
            const traces = [];
            
            // IR (inverted transmittance style)
            traces.push({{
                x: irBroad.x,
                y: irNorm.map(v => 1 - v),
                mode: 'lines',
                line: {{color: '#EF553B', width: 1.5}},
                name: 'IR (Transmittance)',
                hovertemplate: 'IR: %{{x:.0f}} cm⁻¹<extra></extra>'
            }});
            
            // Raman (shifted down)
            const ramanShift = -0.3;
            traces.push({{
                x: ramanBroad.x,
                y: ramanNorm.map(v => v + ramanShift - 1),
                mode: 'lines',
                line: {{color: '#636EFA', width: 1.5}},
                name: 'Raman (Intensity)',
                hovertemplate: 'Raman: %{{x:.0f}} cm⁻¹<extra></extra>'
            }});
            
            // Connector lines
            if (showConnectors) {{
                pairs.forEach((p, i) => {{
                    const irY = 1 - p.irInt;
                    const ramanY = p.ramanInt + ramanShift - 1;
                    const midX = (p.irFreq + p.ramanFreq) / 2;
                    const midY = (irY + ramanY) / 2;
                    
                    traces.push({{
                        x: [p.irFreq, p.irFreq, midX, p.ramanFreq, p.ramanFreq],
                        y: [irY, 0.1, midY, ramanShift - 0.1, ramanY],
                        mode: 'lines',
                        line: {{color: 'gray', width: 1, dash: 'dash'}},
                        showlegend: false,
                        hoverinfo: 'skip'
                    }});
                }});
            }}
            
            const layout = {{
                title: `IR-Raman Correlation: ${{mol.label}}`,
                xaxis: {{
                    title: 'Wavenumber (cm⁻¹)',
                    autorange: invertX ? 'reversed' : true
                }},
                yaxis: {{
                    tickvals: [1, 0.5, 0, ramanShift - 0.5, ramanShift - 1],
                    ticktext: ['0%', '50%', '100% IR', '0.5', '1.0 Raman'],
                    showgrid: false
                }},
                shapes: [{{
                    type: 'line',
                    x0: freqMin, x1: freqMax,
                    y0: ramanShift, y1: ramanShift,
                    line: {{color: 'black', width: 1}}
                }}],
                paper_bgcolor: '{"#2a2a3e" if dark_theme else "#fff"}',
                plot_bgcolor: '{"#2a2a3e" if dark_theme else "#fff"}',
                font: {{color: '{text_primary}'}},
                showlegend: true,
                legend: {{yanchor: 'top', y: 0.99, xanchor: 'left', x: 0.01}},
                margin: {{t: 50, b: 50, l: 60, r: 30}}
            }};
            
            Plotly.newPlot(div, traces, layout, {{responsive: true}});
            
            // Update stats
            document.getElementById('corr-ir-peaks').textContent = irPeaks.length;
            document.getElementById('corr-raman-peaks').textContent = ramanPeaks2.length;
            document.getElementById('corr-paired').textContent = pairs.length;
        }}
        
        // Generate checkboxes with Select All / Clear All
        function generateCheckboxes(containerId, onChangeHandler) {{
            const div = document.getElementById(containerId);
            if (!div) return;
            
            div.innerHTML = '';
            
            // Add action buttons
            const actionsDiv = document.createElement('div');
            actionsDiv.className = 'multiselect-actions-bar';
            actionsDiv.innerHTML = `
                <button type="button" onclick="selectAllCheckboxes('${{containerId}}', true)">Select All</button>
                <button type="button" onclick="selectAllCheckboxes('${{containerId}}', false)">Clear All</button>
            `;
            div.appendChild(actionsDiv);
            
            // Add checkboxes
            const checkboxContainer = document.createElement('div');
            checkboxContainer.style.display = 'flex';
            checkboxContainer.style.flexWrap = 'wrap';
            checkboxContainer.style.gap = '6px';
            checkboxContainer.style.marginTop = '10px';
            
            molecules.forEach((mol, idx) => {{
                const label = document.createElement('label');
                label.className = 'multiselect-item';
                
                const input = document.createElement('input');
                input.type = 'checkbox';
                input.value = idx;
                input.checked = idx < Math.min(10, molecules.length); // Default first 10
                input.onchange = () => onChangeHandler();
                
                label.appendChild(input);
                label.appendChild(document.createTextNode(mol.label));
                checkboxContainer.appendChild(label);
            }});
            
            div.appendChild(checkboxContainer);
        }}
        
        function selectAllCheckboxes(containerId, checked) {{
            const div = document.getElementById(containerId);
            if (!div) return;
            
            div.querySelectorAll('input[type="checkbox"]').forEach(cb => {{
                cb.checked = checked;
            }});
            
            // Trigger re-render based on container
            if (containerId.includes('spectra')) renderSpectra();
            else if (containerId.includes('energy')) renderEnergyChart();
            else if (containerId.includes('orbital')) renderOrbitals(currentMolIndex);
        }}
    </script>
</body>
</html>'''
    
    return html


def prepare_molecules_json(df: pd.DataFrame) -> str:
    """Prepare comprehensive molecule data as JSON for embedding in HTML."""
    molecules = []
    
    def sanitize(val):
        """Convert NaN/None to JSON-safe None."""
        if val is None:
            return None
        if isinstance(val, float) and (pd.isna(val) or np.isnan(val)):
            return None
        if pd.isna(val):
            return None
        return val
    
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        label = f"{mol_id} [{state}]" if state and str(state) != "nan" else str(mol_id)
        
        mol_data = {
            "label": label,
            "mol_id": str(mol_id) if mol_id else "unknown",
            "state": str(state) if state and str(state) != "nan" else None,
            "smiles": str(row.get("smiles")) if row.get("smiles") and str(row.get("smiles")) != "nan" else None,
            "energy": sanitize(row.get("gibbs_Eh")) or sanitize(row.get("single_point_Eh")),
            "single_point": sanitize(row.get("single_point_Eh")),
            "gibbs": sanitize(row.get("gibbs_Eh")),
            "homo": sanitize(row.get("homo_energy")),
            "lumo": sanitize(row.get("lumo_energy")),
            "gap": sanitize(row.get("homo_lumo_gap")),
            "method": str(row.get("method_id")) if row.get("method_id") and str(row.get("method_id")) != "nan" else None,
            "functional": str(row.get("functional")) if row.get("functional") and str(row.get("functional")) != "nan" else None,
            "basis_set": str(row.get("basis_set")) if row.get("basis_set") and str(row.get("basis_set")) != "nan" else None,
            "charge": sanitize(row.get("charge")),
            "multiplicity": sanitize(row.get("multiplicity")),
            "xyz": None,
            "atoms": 0,
            "coords": [],  # Full coordinate table
            "ir": [],
            "raman": [],
            "tddft": [],  # UV-Vis states
            "mulliken": [],  # Mulliken charges
            "orbitals": []  # Orbital energies
        }
        
        # Build XYZ string
        coords = row.get("cart_coords")
        if coords is not None and hasattr(coords, 'empty') and not coords.empty:
            n_atoms = len(coords)
            mol_data["atoms"] = n_atoms
            xyz_lines = [str(n_atoms), mol_id]
            
            for _, atom in coords.iterrows():
                el = str(atom.get("atom", atom.get("element", "C")))
                x = float(atom.get("x", 0))
                y = float(atom.get("y", 0))
                z = float(atom.get("z", 0))
                xyz_lines.append(f"{el} {x:.6f} {y:.6f} {z:.6f}")
            
            mol_data["xyz"] = "\n".join(xyz_lines)
        
        # Add IR data
        ir = row.get("ir")
        if ir is not None and hasattr(ir, 'empty') and not ir.empty:
            ir_data = []
            for _, peak in ir.iterrows():
                freq_col = [c for c in ir.columns if 'freq' in c.lower()][0] if any('freq' in c.lower() for c in ir.columns) else None
                int_col = [c for c in ir.columns if 'intensity' in c.lower() or 'eps' in c.lower()][0] if any('intensity' in c.lower() or 'eps' in c.lower() for c in ir.columns) else None
                
                if freq_col and int_col:
                    ir_data.append({
                        "freq": float(peak[freq_col]),
                        "intensity": float(peak[int_col])
                    })
            mol_data["ir"] = ir_data
        
        # Add Raman data
        raman = row.get("raman")
        if raman is not None and hasattr(raman, 'empty') and not raman.empty:
            raman_data = []
            for _, peak in raman.iterrows():
                freq_col = [c for c in raman.columns if 'freq' in c.lower()][0] if any('freq' in c.lower() for c in raman.columns) else None
                act_col = [c for c in raman.columns if 'activity' in c.lower()][0] if any('activity' in c.lower() for c in raman.columns) else None
                
                if freq_col and act_col:
                    raman_data.append({
                        "freq": float(peak[freq_col]),
                        "activity": float(peak[act_col])
                    })
            mol_data["raman"] = raman_data
        
        # Add TDDFT data (UV-Vis)
        tddft = row.get("tddft_states")
        if tddft is not None and hasattr(tddft, 'empty') and not tddft.empty:
            tddft_data = []
            for _, state in tddft.iterrows():
                state_entry = {"state": int(state.get("State", 0)) if state.get("State") else 0}
                if "Energy_eV" in tddft.columns:
                    state_entry["eV"] = float(state["Energy_eV"])
                if "Wavelength_nm" in tddft.columns:
                    state_entry["nm"] = float(state["Wavelength_nm"])
                if "f" in tddft.columns:
                    state_entry["f"] = float(state["f"])
                tddft_data.append(state_entry)
            mol_data["tddft"] = tddft_data
        
        # Add Mulliken charges
        mulliken = row.get("mulliken")
        if mulliken is not None and hasattr(mulliken, 'empty') and not mulliken.empty:
            mulliken_data = []
            for _, atom in mulliken.iterrows():
                atom_entry = {}
                if "Atom" in mulliken.columns:
                    atom_entry["atom"] = str(atom["Atom"])
                if "Charge" in mulliken.columns:
                    atom_entry["charge"] = float(atom["Charge"])
                if "Spin" in mulliken.columns:
                    atom_entry["spin"] = float(atom["Spin"]) if pd.notna(atom["Spin"]) else 0
                mulliken_data.append(atom_entry)
            mol_data["mulliken"] = mulliken_data
        
        # Add orbitals data - robust extraction
        orbitals = row.get("orbitals")
        orb_data = []
        if orbitals is not None and hasattr(orbitals, 'empty') and not orbitals.empty:
            # Determine column names
            e_col = next((c for c in ["eV", "Energy_eV", "Energy", "energy"] if c in orbitals.columns), None)
            occ_col = next((c for c in ["Occupation", "OCC", "occupancy"] if c in orbitals.columns), None)
            spin_col = next((c for c in ["Spin", "spin"] if c in orbitals.columns), None)
            eh_col = next((c for c in ["Eh", "Energy_Eh"] if c in orbitals.columns), None)
            
            for _, orb in orbitals.iterrows():
                entry = {}
                # Energy
                if e_col:
                    entry["energy"] = float(orb[e_col])
                elif eh_col:
                    entry["energy"] = float(orb[eh_col]) * 27.2114
                else:
                    continue # Skip if no energy
                
                # Occupation
                if occ_col:
                    entry["occ"] = float(orb[occ_col])
                
                # Spin
                if spin_col:
                    entry["spin"] = str(orb[spin_col])
                
                orb_data.append(entry)
        mol_data["orbitals"] = orb_data
        
        # Add full coordinates table
        coords = row.get("cart_coords")
        if coords is not None and hasattr(coords, 'empty') and not coords.empty:
            for _, atom in coords.iterrows():
                coord_entry = {
                    "el": str(atom.get("atom", atom.get("element", "C"))),
                    "x": float(atom.get("x", 0)),
                    "y": float(atom.get("y", 0)),
                    "z": float(atom.get("z", 0))
                }
                mol_data["coords"].append(coord_entry)

        # 5. 2D Image (SVG) extraction
        svg = None
        smiles = row.get("smiles")
        if smiles and str(smiles) != "nan":
            try:
                from rdkit import Chem
                from rdkit.Chem import AllChem
                from rdkit.Chem.Draw import rdMolDraw2D
                import base64
                
                mol = Chem.MolFromSmiles(str(smiles))
                if mol:
                    # Clean up molecule for drawing
                    to_draw = Chem.RemoveHs(mol)
                    AllChem.Compute2DCoords(to_draw)
                    
                    # Create SVG drawer
                    drawer = rdMolDraw2D.MolDraw2DSVG(450, 300)
                    op = drawer.drawOptions()
                    op.addAtomIndices = False
                    op.bondLineWidth = 2
                    
                    drawer.DrawMolecule(to_draw)
                    drawer.FinishDrawing()
                    svg_str = drawer.GetDrawingText()
                    
                    # Encode to base64 to avoid JSON escaping issues with SVG XML
                    svg = base64.b64encode(svg_str.encode('utf-8')).decode('utf-8')
            except Exception:
                pass
        mol_data["svg"] = svg
        
        molecules.append(mol_data)
    
    return json.dumps(molecules, default=str)
