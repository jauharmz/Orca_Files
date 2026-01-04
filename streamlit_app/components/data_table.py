"""
Data Table Component - Display and filter molecule data

Features:
- Interactive DataFrame display
- Column selection
- Export options
"""

import streamlit as st
import pandas as pd
from typing import Optional
import io


def render_data_table(df: pd.DataFrame):
    """
    Render data table with filtering and export.
    
    Args:
        df: DataFrame with molecule data
    """
    
    st.subheader("📊 Data Table")
    
    if df.empty:
        st.warning("No data available")
        return
    
    # Column selector
    all_columns = df.columns.tolist()
    
    # Default display columns
    default_cols = [
        "molecule_id", "optimized_state", "filename", 
        "gibbs_Eh", "single_point_Eh", "homo_energy", "lumo_energy",
        "functional", "basis_set", "calc_class"
    ]
    available_defaults = [c for c in default_cols if c in all_columns]
    
    with st.expander("⚙️ Column Selection", expanded=False):
        selected_cols = st.multiselect(
            "Columns to Display",
            all_columns,
            default=available_defaults if available_defaults else all_columns[:10],
            key="data_table_cols"
        )
    
    if not selected_cols:
        selected_cols = available_defaults if available_defaults else all_columns[:10]
    
    # Display filtered columns only (exclude complex nested data)
    scalar_cols = []
    for col in selected_cols:
        if col in df.columns:
            # Check if column contains DataFrames or complex objects
            sample = df[col].dropna().head(1)
            if len(sample) == 0 or not isinstance(sample.iloc[0], pd.DataFrame):
                scalar_cols.append(col)
    
    if not scalar_cols:
        st.warning("No displayable columns selected")
        return
    
    display_df = df[scalar_cols].copy()
    
    # Column config for better formatting
    column_config = {}
    for col in display_df.columns:
        if any(kw in col.lower() for kw in ["energy", "_eh", "_ev"]):
            column_config[col] = st.column_config.NumberColumn(col, format="%.6f")
        elif any(kw in col.lower() for kw in ["freq", "cm"]):
            column_config[col] = st.column_config.NumberColumn(col, format="%.2f")
    
    # Stats row
    c1, c2, c3 = st.columns(3)
    with c1:
        st.caption(f"📊 {len(display_df)} rows × {len(scalar_cols)} columns")
    with c2:
        if "molecule_id" in display_df.columns:
            st.caption(f"🧬 {display_df['molecule_id'].nunique()} unique molecules")
    with c3:
        if "optimized_state" in display_df.columns:
            states = display_df["optimized_state"].value_counts()
            st.caption("States: " + ", ".join([f"{s}={c}" for s, c in states.items()]))
    
    # Display table
    st.dataframe(
        display_df,
        use_container_width=True,
        height=500,
        column_config=column_config
    )
    
    # Export options
    with st.expander("📤 Export Options"):
        col1, col2, col3 = st.columns(3)
        
        with col1:
            csv = display_df.to_csv(index=False)
            st.download_button(
                "📥 Download CSV",
                csv,
                "orca_data.csv",
                "text/csv",
                key="data_table_csv"
            )
        
        with col2:
            try:
                buffer = io.BytesIO()
                with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
                    display_df.to_excel(writer, index=False, sheet_name='Data')
                buffer.seek(0)
                st.download_button(
                    "📥 Download Excel",
                    buffer,
                    "orca_data.xlsx",
                    "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    key="data_table_xlsx"
                )
            except ImportError:
                st.caption("Install openpyxl for Excel export")
        
        with col3:
            json_str = display_df.to_json(orient='records', indent=2)
            st.download_button(
                "📥 Download JSON",
                json_str,
                "orca_data.json",
                "application/json",
                key="data_table_json"
            )
