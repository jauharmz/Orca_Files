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
    report_name = st.text_input("Report Filename", value="orca_report.html", key="html_filename")
    if not report_name.endswith(".html"):
        report_name += ".html"
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        generate_clicked = st.button("🚀 Generate Report", type="primary")
    
    if generate_clicked:
        try:
            import os
            from pathlib import Path
            
            # Use script directory for reliable path
            script_dir = Path(__file__).parent.parent  # streamlit_app folder
            save_path = script_dir / report_name
            
            st.caption(f"📂 Target save path: `{save_path}`")
            
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
                
                # Store in session state for download button
                st.session_state['report_html'] = html
                st.session_state['report_name'] = report_name
                
                # Save to disk
                with open(save_path, "w", encoding="utf-8") as f:
                    f.write(html)
                
                file_size_mb = len(html) / (1024 * 1024)
                
            st.success(f"✅ Report generated! ({file_size_mb:.1f} MB)")
            
            # Immediate Verification
            if save_path.exists():
                st.info(f"💾 File SAVED to: `{save_path}`")
                st.caption(f"Size: {save_path.stat().st_size / 1024:.1f} KB")
            else:
                st.error(f"❌ File NOT FOUND at expected path!")
            
        except Exception as e:
            st.error(f"Failed to generate report: {e}")
            import traceback
            st.code(traceback.format_exc())

    # Show download button if report exists in session state (always visible after generation)
    if 'report_html' in st.session_state:
        with col2:
            st.download_button(
                label="📥 Download HTML Report",
                data=st.session_state['report_html'],
                file_name=st.session_state.get('report_name', 'orca_report.html'),
                mime="text/html",
                type="primary",
                key="html_download_btn"
            )
        st.markdown("---")
        st.caption("💡 If download button doesn't work, check the saved file path above.")


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
            grid-template-columns: repeat(4, 1fr);
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
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
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
        
        /* Coordinate/Data Tables */
        .data-scroll {{
            max-height: 400px;
            overflow-y: auto;
            border-radius: 8px;
            margin-top: 15px;
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
            .stats-grid {{ grid-template-columns: repeat(2, 1fr); }}
            .mol-info-grid {{ grid-template-columns: 1fr; }}
            .controls {{ flex-direction: column; }}
            #viewer-3d {{ height: 350px; }}
            .chart-container {{ height: 350px; }}
            .tabs {{ flex-wrap: wrap; }}
            .tab {{ flex: 1 1 auto; text-align: center; padding: 10px; }}
            table {{ font-size: 12px; }}
            td, th {{ padding: 8px 6px; }}
            .settings-row {{ flex-direction: column; }}
        }}
        
        @media (max-width: 480px) {{
            .report-header h1 {{ font-size: 1.3em; }}
            .section {{ padding: 15px 10px; }}
            #viewer-3d {{ height: 280px; }}
            .chart-container {{ height: 280px; }}
            .tab {{ padding: 8px 5px; font-size: 12px; }}
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
                <li><a href="#energy-analysis">4. Energy Analysis</a></li>
                <li><a href="#electronic-structure">5. Electronic Structure</a></li>
                <li><a href="#vibrational-analysis">6. Vibrational Analysis</a></li>
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
                    <div style="display: flex; justify-content: center; align-items: center; background: {"#ffffff" if not dark_theme else "#ffffff"}; border-radius: 8px; padding: 20px; min-height: 400px;">
                        <img id="img-2d" src="" alt="2D Structure" style="max-width: 100%; max-height: 400px; display: none;">
                        <div id="no-2d-msg" style="color: #666;">No 2D image available (requires SMILES or RDKit)</div>
                    </div>
                </div>
'''
    
    html += '''            
                <!-- Molecule Info -->
                <div id="mol-info" class="mol-info-grid"></div>
            </div>
        </section>
'''
    
    # Note: Energy Analysis section is now only included once (after Vibrational Analysis) with full controls
    
    if include_orbitals:
        html += '''        
        <!-- 5. Electronic Structure -->
        <section id="electronic-structure" class="section">
            <h2>🔮 5. Electronic Structure</h2>
            
            <p>
                HOMO (Highest Occupied Molecular Orbital) and LUMO (Lowest Unoccupied Molecular Orbital) energy levels 
                provide insight into the electronic properties, reactivity, and optical characteristics of the molecules.
            </p>
            
            <!-- Orbital Settings -->
            <details class="settings-panel">
                <summary>Orbital Visualization Settings</summary>
                <div class="settings-row">
                    <div class="setting-group">
                        <label>Orbitals to Show (+/- HOMO)</label>
                        <input type="range" id="orb-range-slider" min="3" max="20" value="6" onchange="updateOrbitalRange(this.value)">
                        <span class="range-value" id="orb-range-val">6</span>
                    </div>
                </div>
            </details>
            <div id="orbital-chart" class="chart-container"></div>
        </section>
'''
    
    if include_spectra:
        html += '''        
        <!-- 6. Vibrational Analysis -->
        <section id="vibrational-analysis" class="section">
            <h2>📈 6. Vibrational Analysis</h2>
            
            <p>
                Infrared (IR) and Raman spectra provide fingerprints of molecular vibrations. 
                Use the checkboxes below to overlay multiple spectra.
            </p>
            
            <!-- Spectra Settings -->
            <details class="settings-panel">
                <summary>Spectra Visualization Settings</summary>
                <div class="settings-row">
                    <div class="setting-group">
                        <label>Broadening (FWHM)</label>
                        <input type="range" id="fwhm-slider" min="1" max="50" value="10" onchange="updateSpectraBroadening(this.value)">
                        <span class="range-value" id="fwhm-val">10 cm⁻¹</span>
                    </div>
                    <div class="setting-group">
                        <label>Visualization Mode</label>
                        <select id="spectra-mode" onchange="updateSpectraMode(this.value)">
                             <option value="line">Line (Smoothed)</option>
                             <option value="stick">Stick (Discrete)</option>
                        </select>
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
                <div class="tab active" onclick="showSpectraTab('ir')">IR Spectrum</div>
                <div class="tab" onclick="showSpectraTab('raman')">Raman Spectrum</div>
                <div class="tab" onclick="showSpectraTab('uvvis')">UV-Vis (TDDFT)</div>
            </div>
            <div id="spectra-ir" class="tab-content active">
                <div id="ir-chart" class="chart-container"></div>
            </div>
            <div id="spectra-raman" class="tab-content">
                <div id="raman-chart" class="chart-container"></div>
            </div>
            <div id="spectra-uvvis" class="tab-content">
                <div id="uvvis-chart" class="chart-container"></div>
            </div>
        </section>
'''
    
    if include_energy:
        html += '''
        <!-- 4. Energy Analysis -->
        <section id="energy-analysis" class="section">
            <h2>⚡ 4. Energy Analysis</h2>
            
            <details class="settings-panel">
                <summary>Energy Comparison Settings</summary>
                <div class="setting-group" style="width: 100%;">
                    <label><strong>Select Molecules to Compare:</strong></label>
                    <div id="energy-multiselect" class="multiselect-grid">
                        <!-- Checkboxes injected by JS -->
                    </div>
                </div>
            </details>
            
            <div id="energy-chart" class="chart-container"></div>
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
            document.getElementById('mol-view-3d').classList.remove('active');
            document.getElementById('mol-view-2d').classList.remove('active');
            
            // Deactivate all mol view tabs
            document.querySelectorAll('.viewer-container .tab').forEach(t => t.classList.remove('active'));
            
            // Activate clicked tab
            event.target.classList.add('active');
            
            // Show content
            document.getElementById('mol-view-' + mode).classList.add('active');
        }}
        
        function showSpectraTab(type) {{
            document.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
            document.querySelectorAll('.tab-content').forEach(c => c.classList.remove('active'));
            
            event.target.classList.add('active');
            document.getElementById('spectra-' + type).classList.add('active');
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
            if (irDiv && irDiv.offsetParent !== null) {{
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
                    xaxis: {{title: 'Wavenumber (cm⁻¹)', autorange: 'reversed'}},
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
            if (ramanDiv && ramanDiv.offsetParent !== null) {{
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
            if (uvDiv && uvDiv.offsetParent !== null) {{
                 const traces = [];
                 indices.forEach((idx, i) => {{
                    const mol = molecules[idx];
                    if (!mol.tddft || mol.tddft.length === 0) return;
                    
                    const color = colors[i % colors.length];
                    
                    // Simple line logic for UV (Gaussian broadening on energy/nm)
                    // For now, use stick
                    const x = [];
                    const y = [];
                    mol.tddft.forEach(state => {{
                        if(state.nm && state.f) {{
                            x.push(state.nm, state.nm, null);
                            y.push(0, state.f, null);
                        }}
                    }});
                    
                    traces.push({{
                        x: x, y: y,
                        type: 'scatter', mode: 'lines',
                        line: {{color: color, width: 2}},
                        name: mol.label
                    }});
                 }});

                 Plotly.newPlot('uvvis-chart', traces, {{
                    title: 'UV-Vis Spectrum',
                    xaxis: {{title: 'Wavelength (nm)'}},
                    yaxis: {{title: 'Oscillator Strength'}},
                    paper_bgcolor: '{"#2a2a3e" if dark_theme else "#fff"}',
                    plot_bgcolor: '{"#2a2a3e" if dark_theme else "#fff"}',
                    font: {{color: '{text_primary}'}},
                    showlegend: indices.length > 1,
                    margin: {{t: 30, b: 30, l: 50, r: 20}}
                }}, {{responsive: true}});
            }}
        }}
            
            // Raman Chart
            if (mol.raman && mol.raman.length > 0 && document.getElementById('raman-chart')) {{
                let trace = {{}};
                
                if (spectraMode === 'line') {{
                    const freqs = mol.raman.map(p => p.freq);
                    const minFreq = Math.min(...freqs) - 200;
                    const maxFreq = Math.max(...freqs) + 200;
                    // Remap activity to intensity for broadening function
                    const peaks = mol.raman.map(p => ({{freq: p.freq, intensity: p.activity}}));
                    const broadened = gaussianBroadening(peaks, spectraFWHM, Math.max(0, minFreq), maxFreq);
                    
                    trace = {{
                        x: broadened.x, y: broadened.y,
                        type: 'scatter', mode: 'lines',
                        line: {{color: '#ff6b6b', width: 2}},
                        fill: 'tozeroy', fillcolor: 'rgba(255,107,107,0.1)',
                        name: 'Raman Spectrum'
                    }};
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
                        line: {{color: '#ff6b6b', width: 2}},
                        name: 'Peaks'
                    }};
                }}
                
                Plotly.newPlot('raman-chart', [trace], {{
                    title: 'Raman Spectrum',
                    xaxis: {{title: 'Wavenumber (cm⁻¹)', autorange: 'reversed'}},
                    yaxis: {{title: 'Activity'}},
                    paper_bgcolor: '{"#2a2a3e" if dark_theme else "#fff"}',
                    plot_bgcolor: '{"#2a2a3e" if dark_theme else "#fff"}',
                    font: {{color: '{text_primary}'}}
                }}, {{responsive: true}});
            }}
        }}
            
            

        
        function renderEnergyChart() {{
            const indices = getSelectedIndices('energy-multiselect');
            
            const x = [];
            const y = [];
            const markerColors = [];
            const text = [];
            
            indices.forEach(idx => {{
                const mol = molecules[idx];
                x.push(mol.label);
                y.push(mol.energy || 0);
                markerColors.push(idx === currentMolIndex ? '#667eea' : '#a3bffa');
                text.push(`${{mol.energy ? mol.energy.toFixed(5) : 'N/A'}} Eh`);
            }});
            
            const layout = {{
                title: 'Energy Comparison',
                xaxis: {{title: 'Molecule'}},
                yaxis: {{title: 'Energy (Eh)'}},
                paper_bgcolor: '{"#2a2a3e" if dark_theme else "#fff"}',
                plot_bgcolor: '{"#2a2a3e" if dark_theme else "#fff"}',
                font: {{color: '{text_primary}'}},
                margin: {{l: 60, r: 30, t: 50, b: 80}}
            }};
            
            Plotly.newPlot('energy-chart', [{{
                x: x, y: y,
                type: 'bar',
                marker: {{color: markerColors}},
                text: text,
                hoverinfo: 'x+text'
            }}], layout, {{responsive: true}});
        }}
        
        function updateMolInfo(index) {{
            const mol = molecules[index];
            if (!mol) return;
            
            let html = `
                <div class="stat-card"><div class="stat-value">${{mol.atoms || 0}}</div><div class="stat-label">Atoms</div></div>
                <div class="stat-card"><div class="stat-value">${{mol.energy ? mol.energy.toFixed(4) : 'N/A'}}</div><div class="stat-label">Energy (Eh)</div></div>
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
