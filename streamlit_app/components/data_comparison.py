"""
Data Comparison Component - Cross-Molecule Data Comparison

Features:
- Select single or multiple molecules for comparison
- Display all available data side by side
- Horizontal layout: each molecule's data in columns
- Coordinates and all numeric data displayed
- Export to CSV/Excel
"""

import streamlit as st
import pandas as pd
import numpy as np
from typing import List, Dict, Optional
import io


def render_data_comparison(df: pd.DataFrame):
    """Render data comparison table with side-by-side molecule data."""
    
    st.subheader("🔄 Data Comparison")
    st.caption("Compare data across multiple molecules side by side")
    
    if df.empty:
        st.warning("No data available")
        return
    
    # Build molecule options with state labels
    mol_options = []
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        if state and str(state) != "nan":
            label = f"{mol_id} [{state}]"
        else:
            label = mol_id
        mol_options.append({"label": label, "idx": idx, "mol_id": mol_id})
    
    unique_labels = list(dict.fromkeys([m["label"] for m in mol_options]))
    
    if not unique_labels:
        st.warning("No molecules available")
        return
    
    # Molecule selector
    selected = st.multiselect(
        "Select Molecules to Compare",
        unique_labels,
        default=unique_labels[:min(3, len(unique_labels))],
        key="compare_mol_select"
    )
    
    if not selected:
        st.warning("Select at least one molecule to compare")
        return
    
    # Get selected data
    label_to_idx = {m["label"]: m["idx"] for m in mol_options}
    selected_indices = [label_to_idx[l] for l in selected if l in label_to_idx]
    selected_df = df.loc[selected_indices]
    
    st.caption(f"Comparing {len(selected)} molecules")
    
    # Data type selector
    data_types = st.multiselect(
        "Data Types to Compare",
        ["Geometry (Coordinates)", "Energies", "Orbital Energies", "IR Frequencies", "Raman Frequencies", "TDDFT States", "Internal Coordinates"],
        default=["Geometry (Coordinates)", "Energies"],
        key="compare_data_types"
    )
    
    if not data_types:
        st.warning("Select at least one data type")
        return
    
    # Build comparison table based on selected data types
    comparison_tables = {}
    
    for data_type in data_types:
        if data_type == "Geometry (Coordinates)":
            table = build_geometry_comparison(selected_df, selected)
            if table is not None:
                comparison_tables["Geometry"] = table
        
        elif data_type == "Energies":
            table = build_energy_comparison(selected_df, selected)
            if table is not None:
                comparison_tables["Energies"] = table
        
        elif data_type == "Orbital Energies":
            table = build_orbital_comparison(selected_df, selected)
            if table is not None:
                comparison_tables["Orbitals"] = table
        
        elif data_type == "IR Frequencies":
            table = build_ir_comparison(selected_df, selected)
            if table is not None:
                comparison_tables["IR"] = table
        
        elif data_type == "Raman Frequencies":
            table = build_raman_comparison(selected_df, selected)
            if table is not None:
                comparison_tables["Raman"] = table
        
        elif data_type == "TDDFT States":
            table = build_tddft_comparison(selected_df, selected)
            if table is not None:
                comparison_tables["TDDFT"] = table
        
        elif data_type == "Internal Coordinates":
            table = build_internal_coords_comparison(selected_df, selected)
            if table is not None:
                comparison_tables["Internal Coords"] = table
    
    # Display tables
    if not comparison_tables:
        st.info("No data available for the selected data types and molecules")
        return
    
    for name, table_df in comparison_tables.items():
        with st.expander(f"📋 {name}", expanded=True):
            st.dataframe(table_df, use_container_width=True, hide_index=True, height=400)
    
    # Export all tables
    with st.expander("📤 Export All"):
        render_export_combined(comparison_tables)


def build_geometry_comparison(df: pd.DataFrame, labels: List[str]) -> Optional[pd.DataFrame]:
    """Build geometry comparison table with coordinates for each molecule."""
    
    all_rows = []
    max_atoms = 0
    mol_data = {}
    
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        label = f"{mol_id} [{state}]" if state and str(state) != "nan" else mol_id
        
        coords = row.get("optimized_coordinates")
        if coords is not None and hasattr(coords, '__len__') and len(coords) > 0:
            mol_data[label] = coords
            if hasattr(coords, '__len__'):
                max_atoms = max(max_atoms, len(coords))
    
    if not mol_data:
        return None
    
    # Build side-by-side table
    rows = []
    for i in range(max_atoms):
        row_data = {"No": i + 1}
        
        for label in labels:
            if label in mol_data:
                coords = mol_data[label]
                if hasattr(coords, 'iloc') and i < len(coords):  # DataFrame
                    atom_row = coords.iloc[i]
                    if "element" in coords.columns:
                        row_data[f"{label}_Elem"] = atom_row.get("element", "")
                    if "x" in coords.columns:
                        row_data[f"{label}_X"] = f"{atom_row.get('x', 0):.4f}"
                        row_data[f"{label}_Y"] = f"{atom_row.get('y', 0):.4f}"
                        row_data[f"{label}_Z"] = f"{atom_row.get('z', 0):.4f}"
                elif isinstance(coords, list) and i < len(coords):  # List
                    atom = coords[i]
                    if isinstance(atom, dict):
                        row_data[f"{label}_Elem"] = atom.get("element", "")
                        row_data[f"{label}_X"] = f"{atom.get('x', 0):.4f}"
                        row_data[f"{label}_Y"] = f"{atom.get('y', 0):.4f}"
                        row_data[f"{label}_Z"] = f"{atom.get('z', 0):.4f}"
            else:
                row_data[f"{label}_Elem"] = ""
                row_data[f"{label}_X"] = ""
                row_data[f"{label}_Y"] = ""
                row_data[f"{label}_Z"] = ""
        
        rows.append(row_data)
    
    if not rows:
        return None
    
    return pd.DataFrame(rows)


def build_energy_comparison(df: pd.DataFrame, labels: List[str]) -> Optional[pd.DataFrame]:
    """Build energy comparison table."""
    
    rows = []
    energy_cols = ["gibbs_Eh", "single_point_Eh", "enthalpy_Eh", "electronic_energy_Eh", 
                   "zpe_Eh", "thermal_correction_Eh", "homo_energy", "lumo_energy"]
    
    for col in energy_cols:
        row_data = {"Property": col.replace("_", " ").title()}
        
        for idx, row in df.iterrows():
            mol_id = row.get("molecule_id", "unknown")
            state = row.get("optimized_state", "")
            label = f"{mol_id} [{state}]" if state and str(state) != "nan" else mol_id
            
            if label in labels:
                val = row.get(col)
                if val is not None and not pd.isna(val):
                    row_data[label] = f"{val:.6f}"
                else:
                    row_data[label] = "—"
        
        # Only add if at least one molecule has data
        if any(row_data.get(l, "—") != "—" for l in labels):
            rows.append(row_data)
    
    if not rows:
        return None
    
    return pd.DataFrame(rows)


def build_orbital_comparison(df: pd.DataFrame, labels: List[str]) -> Optional[pd.DataFrame]:
    """Build orbital energy comparison table."""
    
    rows = []
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        label = f"{mol_id} [{state}]" if state and str(state) != "nan" else mol_id
        
        if label not in labels:
            continue
        
        orbitals = row.get("orbitals")
        if orbitals is None or (hasattr(orbitals, 'empty') and orbitals.empty):
            continue
        
        homo = row.get("homo_energy")
        
        # Get N closest orbitals to HOMO
        if "energy" in orbitals.columns:
            energies = orbitals["energy"].dropna().sort_values()
        elif "Energy" in orbitals.columns:
            energies = orbitals["Energy"].dropna().sort_values()
        else:
            continue
        
        for i, e in enumerate(energies.values[-10:]):  # Last 10 orbitals
            is_homo = homo is not None and abs(e - homo) < 0.01
            orbital_label = f"{'HOMO' if is_homo else i}"
            
            # Find or create row
            existing = [r for r in rows if r.get("Orbital") == orbital_label]
            if existing:
                existing[0][label] = f"{e:.4f}"
            else:
                row_data = {"Orbital": orbital_label}
                row_data[label] = f"{e:.4f}"
                rows.append(row_data)
    
    if not rows:
        return None
    
    return pd.DataFrame(rows)


def build_ir_comparison(df: pd.DataFrame, labels: List[str]) -> Optional[pd.DataFrame]:
    """Build IR frequency comparison table."""
    
    all_ir_data = {}
    max_rows = 0
    
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        label = f"{mol_id} [{state}]" if state and str(state) != "nan" else mol_id
        
        if label not in labels:
            continue
        
        ir = row.get("ir")
        if ir is not None and hasattr(ir, 'empty') and not ir.empty:
            all_ir_data[label] = ir
            max_rows = max(max_rows, len(ir))
    
    if not all_ir_data:
        return None
    
    rows = []
    for i in range(min(max_rows, 50)):  # Limit to 50 rows
        row_data = {"No": i + 1}
        
        for label in labels:
            if label in all_ir_data:
                ir = all_ir_data[label]
                if i < len(ir):
                    if "freq_cm-1" in ir.columns:
                        row_data[f"{label}_Freq"] = f"{ir.iloc[i]['freq_cm-1']:.1f}"
                    if "intensity_km/mol" in ir.columns:
                        row_data[f"{label}_Int"] = f"{ir.iloc[i]['intensity_km/mol']:.2f}"
                    elif "intensity" in ir.columns:
                        row_data[f"{label}_Int"] = f"{ir.iloc[i]['intensity']:.2f}"
        
        rows.append(row_data)
    
    return pd.DataFrame(rows)


def build_raman_comparison(df: pd.DataFrame, labels: List[str]) -> Optional[pd.DataFrame]:
    """Build Raman frequency comparison table."""
    
    all_raman_data = {}
    max_rows = 0
    
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        label = f"{mol_id} [{state}]" if state and str(state) != "nan" else mol_id
        
        if label not in labels:
            continue
        
        raman = row.get("raman")
        if raman is not None and hasattr(raman, 'empty') and not raman.empty:
            all_raman_data[label] = raman
            max_rows = max(max_rows, len(raman))
    
    if not all_raman_data:
        return None
    
    rows = []
    for i in range(min(max_rows, 50)):
        row_data = {"No": i + 1}
        
        for label in labels:
            if label in all_raman_data:
                raman = all_raman_data[label]
                if i < len(raman):
                    if "freq_cm-1" in raman.columns:
                        row_data[f"{label}_Freq"] = f"{raman.iloc[i]['freq_cm-1']:.1f}"
                    if "activity" in raman.columns:
                        row_data[f"{label}_Act"] = f"{raman.iloc[i]['activity']:.2f}"
        
        rows.append(row_data)
    
    return pd.DataFrame(rows)


def build_tddft_comparison(df: pd.DataFrame, labels: List[str]) -> Optional[pd.DataFrame]:
    """Build TDDFT state comparison table."""
    
    all_tddft_data = {}
    max_states = 0
    
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        label = f"{mol_id} [{state}]" if state and str(state) != "nan" else mol_id
        
        if label not in labels:
            continue
        
        tddft = row.get("tddft_states")
        if tddft is not None and hasattr(tddft, 'empty') and not tddft.empty:
            all_tddft_data[label] = tddft
            if "state" in tddft.columns:
                max_states = max(max_states, tddft["state"].nunique())
    
    if not all_tddft_data:
        return None
    
    rows = []
    for i in range(min(max_states, 10)):
        row_data = {"State": f"S{i+1}"}
        
        for label in labels:
            if label in all_tddft_data:
                tddft = all_tddft_data[label]
                if "state" in tddft.columns:
                    state_data = tddft[tddft["state"] == i + 1]
                    if len(state_data) > 0:
                        s = state_data.iloc[0]
                        if "energy_ev" in tddft.columns:
                            row_data[f"{label}_eV"] = f"{s['energy_ev']:.3f}"
                            row_data[f"{label}_nm"] = f"{1239.84 / s['energy_ev']:.1f}"
                        if "fosc" in tddft.columns:
                            row_data[f"{label}_f"] = f"{s['fosc']:.4f}"
        
        rows.append(row_data)
    
    return pd.DataFrame(rows)


def build_internal_coords_comparison(df: pd.DataFrame, labels: List[str]) -> Optional[pd.DataFrame]:
    """Build internal coordinates comparison table."""
    
    all_data = {}
    all_definitions = set()
    
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        label = f"{mol_id} [{state}]" if state and str(state) != "nan" else mol_id
        
        if label not in labels:
            continue
        
        int_coords = row.get("internal_coords")
        if int_coords is not None and hasattr(int_coords, 'empty') and not int_coords.empty:
            all_data[label] = int_coords
            if "definition" in int_coords.columns:
                all_definitions.update(int_coords["definition"].tolist())
    
    if not all_data or not all_definitions:
        return None
    
    rows = []
    for defn in sorted(all_definitions)[:50]:  # Limit
        row_data = {"Coordinate": defn}
        
        for label in labels:
            if label in all_data:
                coords = all_data[label]
                if "definition" in coords.columns:
                    match = coords[coords["definition"] == defn]
                    if len(match) > 0 and "value" in match.columns:
                        row_data[label] = f"{match['value'].iloc[0]:.4f}"
                    else:
                        row_data[label] = "—"
            else:
                row_data[label] = "—"
        
        rows.append(row_data)
    
    return pd.DataFrame(rows)


def render_export_combined(tables: Dict[str, pd.DataFrame]):
    """Render export options for all comparison tables."""
    
    col1, col2 = st.columns(2)
    
    with col1:
        # CSV export - concatenate all tables
        combined_csv = ""
        for name, table_df in tables.items():
            combined_csv += f"# {name}\n"
            combined_csv += table_df.to_csv(index=False)
            combined_csv += "\n\n"
        
        st.download_button(
            "📥 Download All as CSV",
            combined_csv,
            "data_comparison.csv",
            "text/csv",
            key="compare_csv"
        )
    
    with col2:
        try:
            buffer = io.BytesIO()
            with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
                for name, table_df in tables.items():
                    # Truncate sheet name to 31 chars (Excel limit)
                    sheet_name = name[:31]
                    table_df.to_excel(writer, index=False, sheet_name=sheet_name)
            buffer.seek(0)
            st.download_button(
                "📥 Download All as Excel",
                buffer,
                "data_comparison.xlsx",
                "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                key="compare_xlsx"
            )
        except ImportError:
            st.caption("Install openpyxl for Excel export")
