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
    
    # Initialize session state for HTML
    if "generated_html" not in st.session_state:
        st.session_state.generated_html = None
        st.session_state.html_size = 0
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        if st.button("🚀 Generate HTML Report", type="primary"):
            try:
                with st.spinner("Generating report..."):
                    html = generate_html_report(
                        df, 
                        include_3d=include_3d,
                        include_spectra=include_spectra,
                        include_energy=include_energy,
                        include_orbitals=include_orbitals,
                        dark_theme=dark_theme,
                        compact_mode=compact_mode
                    )
                    st.session_state.generated_html = html
                    st.session_state.html_size = len(html) // 1024
                st.rerun()
                
            except Exception as e:
                st.error(f"Failed to generate report: {e}")
                import traceback
                st.code(traceback.format_exc())
    
    with col2:
        # Show download button if HTML was generated
        if st.session_state.generated_html:
            st.download_button(
                "⬇️ Download HTML Report",
                st.session_state.generated_html,
                "orca_report.html",
                "text/html",
                key="html_download"
            )
            st.success(f"✅ Report ready! ({st.session_state.html_size} KB)")


def generate_html_report(
    df: pd.DataFrame, 
    include_3d: bool = True,
    include_spectra: bool = True,
    include_energy: bool = True,
    include_orbitals: bool = True,
    dark_theme: bool = False,
    compact_mode: bool = False
) -> str:
    """Generate comprehensive standalone HTML report with all visualizations."""
    
    # Theme colors
    if dark_theme:
        bg_primary = "#1e1e2e"
        bg_secondary = "#2a2a3e"
        text_primary = "#ffffff"
        text_secondary = "#a0a0b0"
        accent = "#00ff88"
    else:
        bg_primary = "#f8f9fa"
        bg_secondary = "#ffffff"
        text_primary = "#333333"
        text_secondary = "#666666"
        accent = "#0066cc"
    
    # Prepare molecule data as JSON
    molecules_json = prepare_molecules_json(df)
    
    # Build HTML
    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ORCA Visualization Report</title>
    
    <!-- 3Dmol.js for 3D molecules -->
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    
    <!-- Plotly for charts -->
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    
    <style>
        * {{ box-sizing: border-box; margin: 0; padding: 0; }}
        
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: {bg_primary};
            color: {text_primary};
            line-height: 1.6;
            padding: {"10px" if compact_mode else "20px"};
        }}
        
        .container {{ max-width: 1400px; margin: 0 auto; }}
        
        header {{
            text-align: center;
            padding: {"15px" if compact_mode else "30px"};
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border-radius: 12px;
            margin-bottom: 20px;
        }}
        
        header h1 {{ font-size: {"1.5em" if compact_mode else "2em"}; margin-bottom: 5px; }}
        header p {{ opacity: 0.9; }}
        
        .section {{
            background: {bg_secondary};
            border-radius: 12px;
            padding: {"15px" if compact_mode else "25px"};
            margin-bottom: 20px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        
        .section h2 {{
            color: {accent};
            border-bottom: 2px solid {accent};
            padding-bottom: 10px;
            margin-bottom: 20px;
            font-size: {"1.1em" if compact_mode else "1.3em"};
        }}
        
        .controls {{
            display: flex;
            gap: 15px;
            flex-wrap: wrap;
            margin-bottom: 20px;
            align-items: center;
        }}
        
        select, button {{
            padding: 10px 15px;
            border-radius: 8px;
            border: 1px solid {"#444" if dark_theme else "#ddd"};
            background: {bg_primary};
            color: {text_primary};
            font-size: 14px;
            cursor: pointer;
        }}
        
        button:hover {{ background: {accent}; color: white; }}
        button.active {{ background: {accent}; color: white; }}
        
        #viewer-3d {{
            width: 100%;
            height: 450px;
            border-radius: 8px;
            border: 1px solid {"#444" if dark_theme else "#ddd"};
        }}
        
        .chart-container {{ width: 100%; height: 400px; }}
        
        .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
            font-size: {"12px" if compact_mode else "14px"};
        }}
        
        th, td {{
            padding: {"6px 10px" if compact_mode else "10px 15px"};
            text-align: left;
            border-bottom: 1px solid {"#444" if dark_theme else "#eee"};
        }}
        
        th {{ background: {"#2a2a3e" if dark_theme else "#f0f0f0"}; font-weight: 600; }}
        tr:hover {{ background: {"#333" if dark_theme else "#f8f8f8"}; }}
        
        .metric-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 15px; margin-bottom: 20px; }}
        
        .metric {{
            background: {"#333" if dark_theme else "#f0f5ff"};
            padding: 15px;
            border-radius: 8px;
            text-align: center;
        }}
        
        .metric-value {{ font-size: 1.5em; font-weight: bold; color: {accent}; }}
        .metric-label {{ font-size: 0.85em; color: {text_secondary}; }}
        
        .tabs {{ display: flex; gap: 5px; margin-bottom: 15px; }}
        .tab {{ padding: 8px 16px; border-radius: 6px 6px 0 0; cursor: pointer; }}
        .tab.active {{ background: {accent}; color: white; }}
        
        .tab-content {{ display: none; }}
        .tab-content.active {{ display: block; }}
        
        .tooltip {{
            position: absolute;
            background: rgba(0,0,0,0.9);
            color: white;
            padding: 8px 12px;
            border-radius: 6px;
            font-size: 12px;
            pointer-events: none;
            z-index: 1000;
        }}
        
        footer {{
            text-align: center;
            padding: 20px;
            color: {text_secondary};
            font-size: 12px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>🔬 ORCA Visualization Report</h1>
            <p>Generated by ORCA Visualization Platform</p>
        </header>
        
        <!-- Summary Metrics -->
        <div class="section">
            <div class="metric-grid">
                <div class="metric">
                    <div class="metric-value">{len(df)}</div>
                    <div class="metric-label">Records</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{df["molecule_id"].nunique() if "molecule_id" in df.columns else 0}</div>
                    <div class="metric-label">Molecules</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{df["optimized_state"].nunique() if "optimized_state" in df.columns else 0}</div>
                    <div class="metric-label">States</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{df["gibbs_Eh"].notna().sum() if "gibbs_Eh" in df.columns else 0}</div>
                    <div class="metric-label">With Energy</div>
                </div>
            </div>
        </div>
        
        <!-- Molecule Selector -->
        <div class="section">
            <h2>🧬 Molecule Viewer</h2>
            <div class="controls">
                <label>Select Molecule:</label>
                <select id="molecule-select" onchange="switchMolecule()">
'''
    
    # Add molecule options
    for i, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        label = f"{mol_id} [{state}]" if state and str(state) != "nan" else mol_id
        html += f'                    <option value="{i}">{label}</option>\n'
    
    html += f'''                </select>
'''
    
    if include_3d:
        html += '''                <button onclick="setStyle('stick')">Stick</button>
                <button onclick="setStyle('sphere')">Spacefill</button>
                <button onclick="setStyle('ballstick')" class="active">Ball & Stick</button>
                <button onclick="toggleSpin()">🔄 Spin</button>
'''
    
    html += '''            </div>
'''
    
    if include_3d:
        html += f'''            <div id="viewer-3d" style="background: {"#1e1e2e" if dark_theme else "#f0f0f5"};"></div>
'''
    
    html += '''            
            <!-- Molecule Info -->
            <div id="mol-info" class="metric-grid" style="margin-top: 20px;"></div>
        </div>
'''
    
    if include_spectra:
        html += '''        
        <!-- Spectra Section -->
        <div class="section">
            <h2>📈 Spectra</h2>
            <div class="tabs">
                <div class="tab active" onclick="showSpectraTab('ir')">IR</div>
                <div class="tab" onclick="showSpectraTab('raman')">Raman</div>
                <div class="tab" onclick="showSpectraTab('uv')">UV-Vis</div>
            </div>
            <div id="spectra-ir" class="tab-content active">
                <div id="ir-chart" class="chart-container"></div>
            </div>
            <div id="spectra-raman" class="tab-content">
                <div id="raman-chart" class="chart-container"></div>
            </div>
            <div id="spectra-uv" class="tab-content">
                <div id="uv-chart" class="chart-container"></div>
            </div>
        </div>
'''
    
    if include_energy:
        html += '''        
        <!-- Energy Section -->
        <div class="section">
            <h2>⚡ Energy Comparison</h2>
            <div id="energy-chart" class="chart-container"></div>
        </div>
'''
    
    if include_orbitals:
        html += '''        
        <!-- Orbitals Section -->
        <div class="section">
            <h2>🔮 Molecular Orbitals</h2>
            <div id="orbital-chart" class="chart-container"></div>
        </div>
'''
    
    # Data table
    html += f'''        
        <!-- Data Table -->
        <div class="section">
            <h2>📊 Complete Data</h2>
            <div style="overflow-x: auto;">
                <table>
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
        </div>
        
        <footer>
            Generated by ORCA Visualization Platform | {pd.Timestamp.now().strftime("%Y-%m-%d %H:%M")}
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
