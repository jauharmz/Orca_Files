"""
Data Table Component - Enhanced Display with Comparison Tables

Features:
- Interactive DataFrame display with column selection
- Internal Coordinates comparison (Bonds/Angles/Dihedrals)
- TDDFT State comparison table
- Export options (CSV, Excel, JSON)
"""

import streamlit as st
import pandas as pd
import numpy as np
from typing import Optional, List, Dict
import io


def render_data_table(df: pd.DataFrame):
    """Render data table with filtering, comparison views, and export."""
    
    st.subheader("📊 Data Table")
    
    if df.empty:
        st.warning("No data available")
        return
    
    # Tabs for different views
    data_tabs = st.tabs(["📋 Main Data", "📐 Internal Coordinates", "🔬 TDDFT States"])
    
    with data_tabs[0]:
        render_main_data_table(df)
    with data_tabs[1]:
        render_internal_coords_comparison(df)
    with data_tabs[2]:
        render_tddft_comparison(df)


def render_main_data_table(df: pd.DataFrame):
    """Render main data table with filtering and export."""
    
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
    render_export_options(display_df, "main")


def render_internal_coords_comparison(df: pd.DataFrame):
    """Render internal coordinates comparison table."""
    
    st.markdown("### 📐 Internal Coordinates Comparison")
    st.caption("Compare bonds, angles, and dihedrals across molecules")
    
    # Build molecule options
    mol_options = []
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        if state and str(state) != "nan":
            label = f"{mol_id} [{state}]"
        else:
            label = mol_id
        
        # Check if internal coords available
        int_coords = row.get("internal_coords")
        if int_coords is not None and hasattr(int_coords, 'empty') and not int_coords.empty:
            mol_options.append({"label": label, "idx": idx})
    
    if not mol_options:
        st.info("No internal coordinate data available. This requires geometry optimization output with internal coordinates.")
        return
    
    unique_labels = [m["label"] for m in mol_options]
    
    # Selection
    col1, col2 = st.columns([3, 1])
    with col1:
        selected = st.multiselect(
            "Select Molecules to Compare",
            unique_labels,
            default=unique_labels[:min(4, len(unique_labels))],
            key="int_coords_select"
        )
    with col2:
        coord_filter = st.selectbox(
            "Coordinate Type",
            ["All", "Bonds (B)", "Angles (A)", "Dihedrals (D)"],
            key="coord_filter"
        )
    
    if not selected:
        st.warning("Select molecules to compare")
        return
    
    # Get selected data
    label_to_idx = {m["label"]: m["idx"] for m in mol_options}
    
    # Collect internal coordinates
    all_coords = {}
    for label in selected:
        if label in label_to_idx:
            idx = label_to_idx[label]
            row = df.loc[idx]
            int_coords = row.get("internal_coords")
            if int_coords is not None and hasattr(int_coords, "empty") and not int_coords.empty:
                all_coords[label] = int_coords
    
    if not all_coords:
        st.warning("No internal coordinate data for selected molecules")
        return
    
    # Build comparison table
    comparison_rows = []
    
    # Get all unique coordinate definitions
    all_definitions = set()
    for label, coords in all_coords.items():
        if "definition" in coords.columns:
            all_definitions.update(coords["definition"].tolist())
        elif "type" in coords.columns and "atoms" in coords.columns:
            for _, r in coords.iterrows():
                defn = f"{r.get('type', '')}: {r.get('atoms', '')}"
                all_definitions.add(defn)
    
    # Filter by coordinate type
    filter_map = {
        "Bonds (B)": ["B", "BOND", "STRE"],
        "Angles (A)": ["A", "ANGLE", "BEND"],
        "Dihedrals (D)": ["D", "DIHEDRAL", "TORS"],
        "All": None
    }
    coord_types = filter_map.get(coord_filter)
    
    for defn in sorted(all_definitions):
        row_data = {"Coordinate": defn}
        
        # Check type filter
        if coord_types:
            defn_upper = defn.upper()
            if not any(ct in defn_upper for ct in coord_types):
                continue
        
        for label, coords in all_coords.items():
            value = None
            if "definition" in coords.columns:
                match = coords[coords["definition"] == defn]
                if len(match) > 0 and "value" in match.columns:
                    value = match["value"].iloc[0]
            
            row_data[label] = f"{value:.4f}" if value is not None else "—"
        
        comparison_rows.append(row_data)
    
    if not comparison_rows:
        st.info("No coordinates match the selected filter")
        return
    
    comparison_df = pd.DataFrame(comparison_rows)
    
    st.dataframe(comparison_df, use_container_width=True, height=400, hide_index=True)
    
    # Export
    render_export_options(comparison_df, "internal_coords")


def render_tddft_comparison(df: pd.DataFrame):
    """Render TDDFT state comparison table."""
    
    st.markdown("### 🔬 TDDFT State Comparison")
    st.caption("Compare excited states across molecules")
    
    # Build molecule options
    mol_options = []
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        if state and str(state) != "nan":
            label = f"{mol_id} [{state}]"
        else:
            label = mol_id
        
        # Check if TDDFT data available
        tddft = row.get("tddft_states")
        if tddft is not None and hasattr(tddft, 'empty') and not tddft.empty:
            mol_options.append({"label": label, "idx": idx})
    
    if not mol_options:
        st.info("No TDDFT data available. This requires TDDFT calculation output.")
        return
    
    unique_labels = [m["label"] for m in mol_options]
    
    # Selection
    col1, col2, col3 = st.columns([3, 1, 1])
    with col1:
        selected = st.multiselect(
            "Select Molecules",
            unique_labels,
            default=unique_labels[:min(4, len(unique_labels))],
            key="tddft_select"
        )
    with col2:
        n_states = st.number_input("Number of States", 1, 20, 5, key="tddft_n_states")
    with col3:
        sort_by_fosc = st.checkbox("Sort by f(osc)", True, key="tddft_sort")
    
    if not selected:
        st.warning("Select molecules to compare")
        return
    
    # Get selected data
    label_to_idx = {m["label"]: m["idx"] for m in mol_options}
    
    # Build comparison data
    all_state_data = []
    
    for label in selected:
        if label not in label_to_idx:
            continue
        
        idx = label_to_idx[label]
        row = df.loc[idx]
        tddft = row.get("tddft_states")
        
        if tddft is None or tddft.empty:
            continue
        
        # Get unique states
        if "state" in tddft.columns:
            states_df = tddft.drop_duplicates(subset=["state"]).copy()
        else:
            states_df = tddft.copy()
        
        # Get energy and oscillator strength
        if "energy_ev" in states_df.columns:
            states_df["wavelength_nm"] = 1239.84 / states_df["energy_ev"]
        
        # Sort by oscillator strength if requested
        if sort_by_fosc and "fosc" in states_df.columns:
            states_df = states_df.nlargest(n_states, "fosc")
        elif sort_by_fosc and "weight" in states_df.columns:
            states_df = states_df.nlargest(n_states, "weight")
        else:
            states_df = states_df.head(n_states)
        
        # Sort by state number for display
        if "state" in states_df.columns:
            states_df = states_df.sort_values("state")
        
        # Add to comparison
        for _, state_row in states_df.iterrows():
            state_num = state_row.get("state", "?")
            energy = state_row.get("energy_ev", 0)
            wl = state_row.get("wavelength_nm", 0)
            fosc = state_row.get("fosc", state_row.get("weight", 0))
            
            all_state_data.append({
                "Molecule": label,
                "State": f"S{state_num}",
                "Energy (eV)": f"{energy:.3f}" if energy else "—",
                "λ (nm)": f"{wl:.1f}" if wl else "—",
                "f(osc)": f"{fosc:.4f}" if fosc else "—"
            })
    
    if not all_state_data:
        st.warning("No TDDFT state data for selected molecules")
        return
    
    tddft_df = pd.DataFrame(all_state_data)
    
    # Display as pivot table if multiple molecules
    if len(selected) > 1:
        st.markdown("**Comparison Table**")
        st.dataframe(tddft_df, use_container_width=True, height=400, hide_index=True)
    else:
        st.dataframe(tddft_df, use_container_width=True, height=400, hide_index=True)
    
    # Export
    render_export_options(tddft_df, "tddft")


def render_export_options(df: pd.DataFrame, prefix: str):
    """Render export options for a DataFrame."""
    
    with st.expander("📤 Export Options"):
        col1, col2, col3 = st.columns(3)
        
        with col1:
            # Excel-safe CSV with BOM
            csv = df.to_csv(index=False, encoding='utf-8-sig')
            st.download_button(
                "📥 Download CSV",
                csv,
                f"orca_{prefix}.csv",
                "text/csv",
                key=f"{prefix}_csv"
            )
        
        with col2:
            try:
                buffer = io.BytesIO()
                with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
                    df.to_excel(writer, index=False, sheet_name='Data')
                buffer.seek(0)
                st.download_button(
                    "📥 Download Excel",
                    buffer,
                    f"orca_{prefix}.xlsx",
                    "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    key=f"{prefix}_xlsx"
                )
            except ImportError:
                st.caption("Install openpyxl for Excel export")
        
        with col3:
            json_str = df.to_json(orient='records', indent=2)
            st.download_button(
                "📥 Download JSON",
                json_str,
                f"orca_{prefix}.json",
                "application/json",
                key=f"{prefix}_json"
            )
