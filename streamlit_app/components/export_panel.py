"""
Export Panel Component

Features:
- Molecule + State filter
- Multiple export formats
- Data download
"""

import streamlit as st
import pandas as pd
import json
from io import BytesIO
from typing import List


def render_export_panel(df: pd.DataFrame):
    """Render export panel with molecule+state filter."""
    
    st.subheader("📤 Export Data")
    
    if df.empty:
        st.warning("No data available")
        return
    
    # Build labels
    labels = []
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        label = f"{mol_id} [{state}]" if state else mol_id
        labels.append({"label": label, "idx": idx})
    
    unique_labels = list(dict.fromkeys([l["label"] for l in labels]))
    
    # Filter
    c1, c2 = st.columns([3, 1])
    with c1:
        selected = st.multiselect("Select for Export", unique_labels, 
                                 default=unique_labels, key="export_sel")
    with c2:
        if st.button("✅ All", key="export_all"):
            selected = unique_labels
    
    if not selected:
        st.warning("Select data to export")
        return
    
    # Get data
    sel_indices = [l["idx"] for l in labels if l["label"] in selected]
    sel_df = df.loc[sel_indices]
    
    st.info(f"Selected {len(sel_df)} records for export")
    
    # Export format tabs
    tabs = st.tabs(["📋 CSV", "📊 Excel", "📄 JSON"])
    
    with tabs[0]:
        export_csv(sel_df)
    with tabs[1]:
        export_excel(sel_df)
    with tabs[2]:
        export_json(sel_df)


def export_csv(df: pd.DataFrame):
    """Export to CSV format."""
    
    st.markdown("##### CSV Export")
    
    # Scalar columns only
    scalar_cols = ["molecule_id", "optimized_state", "filename", "smiles", "charge", "multiplicity",
                   "gibbs_Eh", "single_point_Eh", "homo_energy", "lumo_energy", "homo_lumo_gap",
                   "functional", "basis_set", "calc_class"]
    
    available = [c for c in scalar_cols if c in df.columns]
    export_df = df[available].copy()
    
    csv = export_df.to_csv(index=False)
    
    st.download_button(
        "⬇️ Download CSV",
        csv,
        "orca_data.csv",
        "text/csv",
        key="csv_download"
    )
    
    with st.expander("Preview"):
        st.dataframe(export_df.head(10), use_container_width=True)


def export_excel(df: pd.DataFrame):
    """Export to Excel with multiple sheets."""
    
    st.markdown("##### Excel Export (Multiple Sheets)")
    
    output = BytesIO()
    
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        # Main data
        scalar_cols = ["molecule_id", "optimized_state", "filename", "gibbs_Eh", 
                       "single_point_Eh", "homo_energy", "lumo_energy", "calc_class"]
        available = [c for c in scalar_cols if c in df.columns]
        df[available].to_excel(writer, sheet_name="Molecules", index=False)
        
        # Coordinates
        all_coords = []
        for idx, row in df.iterrows():
            mol_id = row.get("molecule_id")
            state = row.get("optimized_state", "")
            coords = row.get("cart_coords")
            if coords is not None and hasattr(coords, 'empty') and not coords.empty:
                coords_copy = coords.copy()
                coords_copy.insert(0, "molecule_id", mol_id)
                coords_copy.insert(1, "optimized_state", state)
                all_coords.append(coords_copy)
        if all_coords:
            pd.concat(all_coords, ignore_index=True).to_excel(writer, sheet_name="Coordinates", index=False)
        
        # Orbitals
        all_orbitals = []
        for idx, row in df.iterrows():
            mol_id = row.get("molecule_id")
            state = row.get("optimized_state", "")
            orbitals = row.get("orbitals")
            if orbitals is not None and hasattr(orbitals, 'empty') and not orbitals.empty:
                orb_copy = orbitals.copy()
                orb_copy.insert(0, "molecule_id", mol_id)
                orb_copy.insert(1, "optimized_state", state)
                all_orbitals.append(orb_copy)
        if all_orbitals:
            pd.concat(all_orbitals, ignore_index=True).to_excel(writer, sheet_name="Orbitals", index=False)
    
    output.seek(0)
    
    st.download_button(
        "⬇️ Download Excel",
        output,
        "orca_data.xlsx",
        "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        key="excel_download"
    )


def export_json(df: pd.DataFrame):
    """Export to JSON format."""
    
    st.markdown("##### JSON Export")
    
    def make_serializable(val):
        if isinstance(val, pd.DataFrame):
            return val.to_dict(orient='records') if not val.empty else None
        if pd.isna(val):
            return None
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
    
    with st.expander("Preview"):
        st.json(records[0] if records else {})
