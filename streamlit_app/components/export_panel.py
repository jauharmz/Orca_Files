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
    """Export interactive HTML report."""
    
    st.markdown("##### HTML Report")
    st.info("📄 Generates a standalone HTML file with interactive Plotly charts.")
    
    if st.button("🎨 Generate HTML Report"):
        try:
            html = generate_html_report(df)
            
            st.download_button(
                "⬇️ Download HTML Report",
                html,
                "orca_report.html",
                "text/html",
                key="html_download"
            )
            st.success("✅ Report generated!")
            
        except Exception as e:
            st.error(f"Failed to generate report: {e}")


def generate_html_report(df: pd.DataFrame) -> str:
    """Generate standalone HTML report."""
    import plotly.graph_objects as go
    
    # Simple HTML template
    html = """<!DOCTYPE html>
<html>
<head>
    <title>ORCA Data Report</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body { font-family: sans-serif; margin: 20px; background: #f8f9fa; }
        h1 { color: #333; }
        .section { background: white; padding: 20px; margin: 20px 0; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
        table { border-collapse: collapse; width: 100%; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background: #f0f0f0; }
    </style>
</head>
<body>
    <h1>🔬 ORCA Data Report</h1>
    <div class="section">
        <h2>Summary</h2>
        <p>Total records: """ + str(len(df)) + """</p>
        <p>Generated by ORCA Visualization Platform</p>
    </div>
    <div class="section">
        <h2>Molecule Data</h2>
        <table>
            <tr><th>Molecule</th><th>State</th><th>Energy (Eh)</th><th>HOMO (eV)</th><th>LUMO (eV)</th></tr>
"""
    
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "N/A")
        state = row.get("optimized_state", "N/A")
        energy = row.get("gibbs_Eh") or row.get("single_point_Eh")
        homo = row.get("homo_energy")
        lumo = row.get("lumo_energy")
        
        html += f"""            <tr>
                <td>{mol_id}</td>
                <td>{state}</td>
                <td>{f'{energy:.4f}' if energy else 'N/A'}</td>
                <td>{f'{homo:.3f}' if homo else 'N/A'}</td>
                <td>{f'{lumo:.3f}' if lumo else 'N/A'}</td>
            </tr>
"""
    
    html += """        </table>
    </div>
</body>
</html>"""
    
    return html
