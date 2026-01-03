"""
Data Table Component - Horizontal multi-molecule layout with export

Creates adaptive tables that show data from multiple molecules side by side.
"""

import streamlit as st
import pandas as pd
from typing import Dict, List, Optional
import io


def create_horizontal_table(
    data_dict: Dict[str, pd.DataFrame],
    coordinate_col: str,
    data_cols: List[str],
    title: str = "Data Table"
) -> pd.DataFrame:
    """
    Create horizontal multi-molecule data table.
    
    Args:
        data_dict: {molecule_id: dataframe}
        coordinate_col: Column to use as coordinate (e.g., "freq_cm-1")
        data_cols: Columns to include as data
        title: Table title
        
    Returns:
        Combined DataFrame with horizontal layout
    """
    if not data_dict:
        return pd.DataFrame()
    
    # Get all unique coordinate values across all molecules
    all_coords = set()
    for mol_id, df in data_dict.items():
        if df is not None and not df.empty and coordinate_col in df.columns:
            all_coords.update(df[coordinate_col].tolist())
    
    all_coords = sorted(all_coords)
    
    if not all_coords:
        return pd.DataFrame()
    
    # Build result DataFrame
    result = pd.DataFrame({"No": range(1, len(all_coords) + 1)})
    
    for mol_id, df in data_dict.items():
        if df is None or df.empty:
            continue
        
        # Add coordinate column for this molecule
        result[f"{mol_id}_{coordinate_col}"] = all_coords
        
        # Add data columns for this molecule
        for col in data_cols:
            if col in df.columns:
                # Map values to coordinates
                coord_to_val = dict(zip(df[coordinate_col], df[col]))
                result[f"{mol_id}_{col}"] = [coord_to_val.get(c, None) for c in all_coords]
    
    return result


def render_data_table(
    data_dict: Dict[str, pd.DataFrame],
    coordinate_col: str,
    data_cols: List[str],
    title: str = "Data Table",
    show_export: bool = True
):
    """
    Render interactive data table with export options.
    
    Args:
        data_dict: {molecule_id: dataframe}
        coordinate_col: Column to use as coordinate
        data_cols: Columns to include as data
        title: Table title
        show_export: Whether to show export buttons
    """
    if not data_dict:
        st.warning("No data to display")
        return
    
    # Create horizontal table
    table_df = create_horizontal_table(data_dict, coordinate_col, data_cols, title)
    
    if table_df.empty:
        st.warning("No data available for selected molecules")
        return
    
    st.subheader(f"📋 {title}")
    
    # Column config for better display
    column_config = {}
    for col in table_df.columns:
        if col == "No":
            column_config[col] = st.column_config.NumberColumn(col, format="%d")
        elif any(keyword in col.lower() for keyword in ["energy", "eh", "ev"]):
            column_config[col] = st.column_config.NumberColumn(col, format="%.6f")
        elif any(keyword in col.lower() for keyword in ["freq", "cm"]):
            column_config[col] = st.column_config.NumberColumn(col, format="%.2f")
        elif any(keyword in col.lower() for keyword in ["intensity", "activity"]):
            column_config[col] = st.column_config.NumberColumn(col, format="%.4f")
    
    # Display table
    st.dataframe(
        table_df,
        use_container_width=True,
        height=400,
        column_config=column_config
    )
    
    # Export options
    if show_export:
        col1, col2, col3 = st.columns(3)
        
        with col1:
            csv = table_df.to_csv(index=False)
            st.download_button(
                "📥 Download CSV",
                csv,
                f"{title.lower().replace(' ', '_')}.csv",
                "text/csv",
                key=f"csv_{title}"
            )
        
        with col2:
            # Excel export
            buffer = io.BytesIO()
            with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
                table_df.to_excel(writer, index=False, sheet_name='Data')
            buffer.seek(0)
            st.download_button(
                "📥 Download Excel",
                buffer,
                f"{title.lower().replace(' ', '_')}.xlsx",
                "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                key=f"xlsx_{title}"
            )
        
        with col3:
            st.write(f"Rows: {len(table_df)} | Cols: {len(table_df.columns)}")
    
    return table_df


def render_simple_table(
    df: pd.DataFrame,
    title: str = "Data",
    show_export: bool = True,
    height: int = 400
):
    """Render a simple DataFrame with export options."""
    
    st.subheader(f"📋 {title}")
    
    st.dataframe(df, use_container_width=True, height=height)
    
    if show_export and not df.empty:
        col1, col2 = st.columns(2)
        
        with col1:
            csv = df.to_csv(index=False)
            st.download_button(
                "📥 CSV",
                csv,
                f"{title.lower().replace(' ', '_')}.csv",
                "text/csv",
                key=f"simple_csv_{title}"
            )
        
        with col2:
            st.write(f"Rows: {len(df)}")
