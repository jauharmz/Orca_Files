"""
Export Panel Component - Export data in various formats

Features:
- CSV export (all data or filtered)
- Excel export with multiple sheets
- JSON export
- HTML report generation
- Filtered data export
"""

import streamlit as st
import pandas as pd
from typing import Optional
import io
import json


def render_export_panel(df: pd.DataFrame):
    """Render export panel."""
    
    st.header("📤 Export Data")
    
    # Data summary
    st.info(f"Data contains {len(df)} molecules")
    
    # Filter options
    st.subheader("🔍 Export Selection")
    
    col1, col2 = st.columns(2)
    
    with col1:
        export_type = st.radio(
            "Export Type",
            ["All Data", "Filtered Data", "Selected Molecules"],
            key="export_type"
        )
    
    with col2:
        if export_type == "Selected Molecules":
            mol_ids = df["molecule_id"].dropna().unique().tolist()
            selected = st.multiselect("Select Molecules", mol_ids, key="export_mol_select")
            export_df = df[df["molecule_id"].isin(selected)] if selected else df
        elif export_type == "Filtered Data":
            # Use global filter from session state
            selected_mols = st.session_state.get("selected_molecules", [])
            if selected_mols:
                export_df = df[df["molecule_id"].isin(selected_mols)]
            else:
                export_df = df
        else:
            export_df = df
    
    st.write(f"Exporting: **{len(export_df)} molecules**")
    
    st.divider()
    
    # Export formats
    st.subheader("📁 Export Formats")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown("### CSV")
        
        # Scalar columns only
        scalar_cols = get_scalar_columns(export_df)
        
        csv_data = export_df[scalar_cols].to_csv(index=False)
        st.download_button(
            "📥 Download CSV",
            csv_data,
            "orca_data.csv",
            "text/csv",
            key="export_csv"
        )
        st.caption("Scalar columns only")
    
    with col2:
        st.markdown("### Excel")
        
        buffer = create_excel_export(export_df)
        st.download_button(
            "📥 Download Excel",
            buffer,
            "orca_data.xlsx",
            "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            key="export_xlsx"
        )
        st.caption("Multiple sheets")
    
    with col3:
        st.markdown("### JSON")
        
        json_data = create_json_export(export_df)
        st.download_button(
            "📥 Download JSON",
            json_data,
            "orca_data.json",
            "application/json",
            key="export_json"
        )
        st.caption("Full nested data")
    
    st.divider()
    
    # Advanced exports
    st.subheader("📊 Advanced Exports")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("**Individual Data Types**")
        
        data_type = st.selectbox(
            "Data Type",
            ["IR Spectra", "Raman Spectra", "TDDFT States", "Coordinates", "Orbitals", "Mulliken"],
            key="export_data_type"
        )
        
        if st.button("📥 Export Selected Type", key="export_type_btn"):
            export_data_type(export_df, data_type)
    
    with col2:
        st.markdown("**Horizontal Table Format**")
        st.caption("Multiple molecules side-by-side")
        
        if st.button("📥 Export Horizontal Table", key="export_horiz_btn"):
            export_horizontal_table(export_df)


def get_scalar_columns(df: pd.DataFrame) -> list:
    """Get scalar (non-nested) columns."""
    scalar_cols = []
    for col in df.columns:
        sample = df[col].dropna().head(1)
        if len(sample) == 0:
            continue
        val = sample.iloc[0]
        if not isinstance(val, (pd.DataFrame, list, dict)):
            scalar_cols.append(col)
    return scalar_cols


def create_excel_export(df: pd.DataFrame) -> bytes:
    """Create Excel file with multiple sheets."""
    buffer = io.BytesIO()
    
    with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
        # Main sheet - scalar data
        scalar_cols = get_scalar_columns(df)
        df[scalar_cols].to_excel(writer, sheet_name='Summary', index=False)
        
        # IR data sheet
        ir_rows = []
        for _, row in df.iterrows():
            ir = row.get("ir")
            if ir is not None and hasattr(ir, 'empty') and not ir.empty:
                ir_copy = ir.copy()
                ir_copy.insert(0, 'molecule_id', row['molecule_id'])
                ir_rows.append(ir_copy)
        if ir_rows:
            ir_df = pd.concat(ir_rows, ignore_index=True)
            ir_df.to_excel(writer, sheet_name='IR_Spectra', index=False)
        
        # Raman data sheet
        raman_rows = []
        for _, row in df.iterrows():
            raman = row.get("raman")
            if raman is not None and hasattr(raman, 'empty') and not raman.empty:
                raman_copy = raman.copy()
                raman_copy.insert(0, 'molecule_id', row['molecule_id'])
                raman_rows.append(raman_copy)
        if raman_rows:
            raman_df = pd.concat(raman_rows, ignore_index=True)
            raman_df.to_excel(writer, sheet_name='Raman_Spectra', index=False)
        
        # TDDFT data sheet
        tddft_rows = []
        for _, row in df.iterrows():
            tddft = row.get("tddft_states")
            if tddft is not None and hasattr(tddft, 'empty') and not tddft.empty:
                tddft_copy = tddft.copy()
                tddft_copy.insert(0, 'molecule_id', row['molecule_id'])
                tddft_rows.append(tddft_copy)
        if tddft_rows:
            tddft_df = pd.concat(tddft_rows, ignore_index=True)
            tddft_df.to_excel(writer, sheet_name='TDDFT_States', index=False)
    
    buffer.seek(0)
    return buffer.getvalue()


def create_json_export(df: pd.DataFrame) -> str:
    """Create JSON export with nested data."""
    
    def serialize_value(val):
        if isinstance(val, pd.DataFrame):
            return val.to_dict(orient='records')
        elif isinstance(val, (pd.Series,)):
            return val.tolist()
        elif pd.isna(val):
            return None
        else:
            return val
    
    records = []
    for _, row in df.iterrows():
        record = {}
        for col in df.columns:
            record[col] = serialize_value(row[col])
        records.append(record)
    
    return json.dumps(records, indent=2, default=str)


def export_data_type(df: pd.DataFrame, data_type: str):
    """Export specific data type."""
    
    type_mapping = {
        "IR Spectra": "ir",
        "Raman Spectra": "raman",
        "TDDFT States": "tddft_states",
        "Coordinates": "cart_coords",
        "Orbitals": "orbitals",
        "Mulliken": "mulliken"
    }
    
    col_name = type_mapping.get(data_type)
    if not col_name:
        st.error("Unknown data type")
        return
    
    rows = []
    for _, row in df.iterrows():
        data = row.get(col_name)
        if data is not None and hasattr(data, 'empty') and not data.empty:
            data_copy = data.copy()
            data_copy.insert(0, 'molecule_id', row['molecule_id'])
            rows.append(data_copy)
    
    if not rows:
        st.warning(f"No {data_type} data available")
        return
    
    combined_df = pd.concat(rows, ignore_index=True)
    
    csv_data = combined_df.to_csv(index=False)
    st.download_button(
        f"📥 Download {data_type}",
        csv_data,
        f"{col_name}.csv",
        "text/csv",
        key=f"download_{col_name}"
    )
    
    st.success(f"Ready to download {len(combined_df)} rows of {data_type}")


def export_horizontal_table(df: pd.DataFrame):
    """Export horizontal multi-molecule table."""
    
    # Select data type
    st.info("Creating horizontal table with IR data...")
    
    ir_data = {}
    for _, row in df.iterrows():
        ir = row.get("ir")
        if ir is not None and hasattr(ir, 'empty') and not ir.empty:
            ir_data[row['molecule_id']] = ir
    
    if not ir_data:
        st.warning("No IR data available")
        return
    
    from components.data_table import create_horizontal_table
    
    table = create_horizontal_table(
        ir_data,
        coordinate_col="freq_cm-1",
        data_cols=["intensity_km/mol"]
    )
    
    if table.empty:
        st.warning("Could not create horizontal table")
        return
    
    csv_data = table.to_csv(index=False)
    st.download_button(
        "📥 Download Horizontal Table",
        csv_data,
        "horizontal_ir_data.csv",
        "text/csv",
        key="download_horiz"
    )
    
    st.success(f"Created table with {len(table)} rows, {len(table.columns)} columns")
