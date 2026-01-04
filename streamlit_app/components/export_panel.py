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
        if st.button("🚀 Generate & Save Report", type="primary"):
            try:
                with st.spinner("Generating comprehensive report..."):
                    html = generate_html_report(
                        df, 
                        include_3d=include_3d,
                        include_spectra=include_spectra,
                        include_energy=include_energy,
                        include_orbitals=include_orbitals,
                        dark_theme=dark_theme,
                        compact_mode=compact_mode
                    )
                    
                    # Save to disk
                    filepath = f"./{report_name}"
                    with open(filepath, "w", encoding="utf-8") as f:
                        f.write(html)
                    
                    file_size = len(html) // 1024
                    
                st.success(f"✅ Report saved to: **{filepath}** ({file_size} KB)")
                st.info(f"📂 Open the file in your browser: `{filepath}`")
                
            except Exception as e:
                st.error(f"Failed to generate report: {e}")
                import traceback
                st.code(traceback.format_exc())


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
        }}
        
        #viewer-3d {{
            width: 100%;
            height: 500px;
            border-radius: 8px;
            border: 1px solid rgba(0,0,0,0.1);
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
            .container {{ padding: 20px 15px; }}
            .section {{ padding: 25px; }}
            .report-header {{ padding: 40px 20px; }}
            .report-header h1 {{ font-size: 1.8em; }}
            .stats-grid {{ grid-template-columns: repeat(2, 1fr); }}
            .mol-info-grid {{ grid-template-columns: repeat(2, 1fr); }}
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
                <div class="controls">
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
'''
    
    if include_3d:
        html += '''                    <button onclick="setStyle('stick')">Stick</button>
                    <button onclick="setStyle('sphere')">Spacefill</button>
                    <button onclick="setStyle('ballstick')" class="active">Ball & Stick</button>
                    <button onclick="toggleSpin()">🔄 Spin</button>
'''
    
    html += '''                </div>
'''
    
    if include_3d:
        html += f'''                <div id="viewer-3d" style="background: {"#1e1e2e" if dark_theme else "#f0f0f5"};"></div>
'''
    
    html += '''            
                <!-- Molecule Info -->
                <div id="mol-info" class="mol-info-grid"></div>
            </div>
        </section>
'''
    
    if include_energy:
        html += '''        
        <!-- 4. Energy Analysis -->
        <section id="energy-analysis" class="section">
            <h2>⚡ 4. Energy Analysis</h2>
            
            <p>
                Comparative energy analysis of all molecular systems. The chart below shows the Gibbs free energy 
                (or single-point energy where thermochemistry is not available) for each calculation.
            </p>
            
            <div id="energy-chart" class="chart-container"></div>
        </section>
'''
    
    if include_orbitals:
        html += '''        
        <!-- 5. Electronic Structure -->
        <section id="electronic-structure" class="section">
            <h2>🔮 5. Electronic Structure</h2>
            
            <p>
                HOMO (Highest Occupied Molecular Orbital) and LUMO (Lowest Unoccupied Molecular Orbital) energy levels 
                provide insight into the electronic properties, reactivity, and optical characteristics of the molecules.
            </p>
            
            <div id="orbital-chart" class="chart-container"></div>
        </section>
'''
    
    if include_spectra:
        html += '''        
        <!-- 6. Vibrational Analysis -->
        <section id="vibrational-analysis" class="section">
            <h2>📈 6. Vibrational Analysis</h2>
            
            <p>
                Infrared (IR) and Raman spectra provide fingerprints of molecular vibrations. Select a molecule above 
                to view its vibrational spectra.
            </p>
            
            <div class="tabs">
                <div class="tab active" onclick="showSpectraTab('ir')">IR Spectrum</div>
                <div class="tab" onclick="showSpectraTab('raman')">Raman Spectrum</div>
            </div>
            <div id="spectra-ir" class="tab-content active">
                <div id="ir-chart" class="chart-container"></div>
            </div>
            <div id="spectra-raman" class="tab-content">
                <div id="raman-chart" class="chart-container"></div>
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
        
        smiles_display = str(smiles)[:30] + "..." if smiles and len(str(smiles)) > 30 else smiles
        
        html += f'''                        <tr>
                            <td><strong>{mol_id}</strong></td>
                            <td>{state}</td>
                            <td title="{smiles}">{smiles_display if smiles and str(smiles) != "nan" else "N/A"}</td>
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
            if (document.getElementById('viewer-3d')) {{
                initViewer();
            }}
            updateMolInfo(0);
            if (document.getElementById('energy-chart')) {{
                renderEnergyChart();
            }}
            if (document.getElementById('orbital-chart')) {{
                renderOrbitalChart();
            }}
            renderSpectra(0);
        }});
        
        function initViewer() {{
            viewer = $3Dmol.createViewer('viewer-3d', {{
                backgroundColor: '{"#1e1e2e" if dark_theme else "#f0f0f5"}'
            }});
            loadMolecule(0);
        }}
        
        function loadMolecule(index) {{
            if (!viewer || !molecules[index]) return;
            currentMolIndex = index;
            
            const mol = molecules[index];
            viewer.removeAllModels();
            
            if (mol.xyz) {{
                viewer.addModel(mol.xyz, 'xyz');
                applyStyle(currentStyle);
                viewer.zoomTo();
                viewer.render();
            }}
            
            updateMolInfo(index);
            renderSpectra(index);
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
            
            if (style === 'stick') {{
                viewer.setStyle({{}}, {{stick: {{radius: 0.15}}}});
            }} else if (style === 'sphere') {{
                viewer.setStyle({{}}, {{sphere: {{scale: 0.8}}}});
            }} else {{
                viewer.setStyle({{}}, {{stick: {{radius: 0.12}}, sphere: {{scale: 0.25}}}});
            }}
            viewer.render();
        }}
        
        function toggleSpin() {{
            spinning = !spinning;
            if (spinning) {{
                spinInterval = setInterval(() => {{
                    viewer.rotate(1, {{x: 0, y: 1, z: 0}});
                    viewer.render();
                }}, 40);
            }} else {{
                clearInterval(spinInterval);
            }}
        }}
        
        function updateMolInfo(index) {{
            const mol = molecules[index];
            if (!mol) return;
            
            const infoDiv = document.getElementById('mol-info');
            infoDiv.innerHTML = `
                <div class="metric"><div class="metric-value">${{mol.atoms || 'N/A'}}</div><div class="metric-label">Atoms</div></div>
                <div class="metric"><div class="metric-value">${{mol.energy ? mol.energy.toFixed(4) : 'N/A'}}</div><div class="metric-label">Energy (Eh)</div></div>
                <div class="metric"><div class="metric-value">${{mol.homo ? mol.homo.toFixed(2) : 'N/A'}}</div><div class="metric-label">HOMO (eV)</div></div>
                <div class="metric"><div class="metric-value">${{mol.lumo ? mol.lumo.toFixed(2) : 'N/A'}}</div><div class="metric-label">LUMO (eV)</div></div>
            `;
        }}
        
        function showSpectraTab(type) {{
            document.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
            document.querySelectorAll('.tab-content').forEach(c => c.classList.remove('active'));
            
            event.target.classList.add('active');
            document.getElementById('spectra-' + type).classList.add('active');
        }}
        
        function renderSpectra(index) {{
            const mol = molecules[index];
            if (!mol) return;
            
            // IR Chart
            if (mol.ir && mol.ir.length > 0 && document.getElementById('ir-chart')) {{
                const x = mol.ir.map(p => p.freq);
                const y = mol.ir.map(p => p.intensity);
                
                Plotly.newPlot('ir-chart', [{{
                    x: x, y: y,
                    type: 'scatter', mode: 'lines',
                    line: {{color: '#00ff88', width: 1.5}},
                    fill: 'tozeroy', fillcolor: 'rgba(0,255,136,0.1)'
                }}], {{
                    title: 'IR Spectrum',
                    xaxis: {{title: 'Wavenumber (cm⁻¹)', autorange: 'reversed'}},
                    yaxis: {{title: 'Intensity'}},
                    paper_bgcolor: '{"#2a2a3e" if dark_theme else "#fff"}',
                    plot_bgcolor: '{"#2a2a3e" if dark_theme else "#fff"}',
                    font: {{color: '{text_primary}'}}
                }}, {{responsive: true}});
            }}
            
            // Raman Chart
            if (mol.raman && mol.raman.length > 0 && document.getElementById('raman-chart')) {{
                const x = mol.raman.map(p => p.freq);
                const y = mol.raman.map(p => p.activity);
                
                Plotly.newPlot('raman-chart', [{{
                    x: x, y: y,
                    type: 'scatter', mode: 'lines',
                    line: {{color: '#ff6b6b', width: 1.5}},
                    fill: 'tozeroy', fillcolor: 'rgba(255,107,107,0.1)'
                }}], {{
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
            const labels = molecules.map(m => m.label);
            const energies = molecules.map(m => m.energy || 0);
            
            Plotly.newPlot('energy-chart', [{{
                x: labels,
                y: energies,
                type: 'bar',
                marker: {{color: '#667eea'}}
            }}], {{
                title: 'Energy Comparison',
                xaxis: {{title: 'Molecule'}},
                yaxis: {{title: 'Energy (Eh)'}},
                paper_bgcolor: '{"#2a2a3e" if dark_theme else "#fff"}',
                plot_bgcolor: '{"#2a2a3e" if dark_theme else "#fff"}',
                font: {{color: '{text_primary}'}}
            }}, {{responsive: true}});
        }}
        
        function renderOrbitalChart() {{
            const labels = molecules.map(m => m.label);
            const homos = molecules.map(m => m.homo || 0);
            const lumos = molecules.map(m => m.lumo || 0);
            
            Plotly.newPlot('orbital-chart', [
                {{x: labels, y: homos, name: 'HOMO', type: 'bar', marker: {{color: '#ff6b6b'}}}},
                {{x: labels, y: lumos, name: 'LUMO', type: 'bar', marker: {{color: '#4ecdc4'}}}}
            ], {{
                title: 'HOMO/LUMO Energy Levels',
                xaxis: {{title: 'Molecule'}},
                yaxis: {{title: 'Energy (eV)'}},
                barmode: 'group',
                paper_bgcolor: '{"#2a2a3e" if dark_theme else "#fff"}',
                plot_bgcolor: '{"#2a2a3e" if dark_theme else "#fff"}',
                font: {{color: '{text_primary}'}}
            }}, {{responsive: true}});
        }}
    </script>
</body>
</html>'''
    
    return html


def prepare_molecules_json(df: pd.DataFrame) -> str:
    """Prepare molecule data as JSON for embedding in HTML."""
    molecules = []
    
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        label = f"{mol_id} [{state}]" if state and str(state) != "nan" else mol_id
        
        mol_data = {
            "label": label,
            "mol_id": mol_id,
            "state": state if state and str(state) != "nan" else None,
            "smiles": row.get("smiles") if row.get("smiles") and str(row.get("smiles")) != "nan" else None,
            "energy": row.get("gibbs_Eh") or row.get("single_point_Eh"),
            "homo": row.get("homo_energy"),
            "lumo": row.get("lumo_energy"),
            "gap": row.get("homo_lumo_gap"),
            "method": row.get("method_id"),
            "xyz": None,
            "atoms": 0,
            "ir": [],
            "raman": []
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
            
            mol_data["xyz"] = "\\n".join(xyz_lines)
        
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
        
        molecules.append(mol_data)
    
    return json.dumps(molecules, default=str)
