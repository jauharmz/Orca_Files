"""
Data Comparison Component - Horizontal Per-Molecule Layout

Features:
- Select single or multiple molecules for comparison
- Horizontal layout: each molecule's data in separate column groups
- Format: |No| coord_mol1 | data1_mol1 | data2_mol1 | coord_mol2 | data1_mol2 | data2_mol2 | ...
- Export to CSV/Excel
"""

import streamlit as st
import pandas as pd
import numpy as np
from typing import List, Dict, Optional
import io


def render_data_comparison(df: pd.DataFrame):
    """Render data comparison table with horizontal per-molecule layout."""
    
    st.subheader("🔄 Data Comparison")
    st.caption("Compare data across multiple molecules - each molecule in its own column group")
    
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
    
    st.caption(f"Comparing {len(selected)} molecules")
    
    # Data type selector
    data_type = st.selectbox(
        "Select Data Type to Compare",
        ["Geometry (Coordinates)", "IR Frequencies", "Raman Frequencies", "TDDFT States", "Internal Coordinates", "Energies"],
        key="compare_data_type"
    )
    
    # Build and display the comparison table
    if data_type == "Geometry (Coordinates)":
        table_df = build_horizontal_geometry_table(df, selected, label_to_idx)
    elif data_type == "IR Frequencies":
        table_df = build_horizontal_ir_table(df, selected, label_to_idx)
    elif data_type == "Raman Frequencies":
        table_df = build_horizontal_raman_table(df, selected, label_to_idx)
    elif data_type == "TDDFT States":
        table_df = build_horizontal_tddft_table(df, selected, label_to_idx)
    elif data_type == "Internal Coordinates":
        table_df = build_horizontal_internal_coords_table(df, selected, label_to_idx)
    elif data_type == "Energies":
        table_df = build_horizontal_energy_table(df, selected, label_to_idx)
    else:
        table_df = None
    
    if table_df is None or table_df.empty:
        st.info(f"No {data_type} data available for selected molecules")
        return
    
    # Display table
    st.dataframe(table_df, use_container_width=True, height=500, hide_index=True)
    
    # Export
    render_export_options(table_df, data_type.lower().replace(" ", "_").replace("(", "").replace(")", ""))


def build_horizontal_geometry_table(df: pd.DataFrame, labels: List[str], label_to_idx: Dict) -> Optional[pd.DataFrame]:
    """Build horizontal geometry table: |No| Element_mol1 | X_mol1 | Y_mol1 | Z_mol1 | Element_mol2 | ... |"""
    
    mol_coords = {}
    max_atoms = 0
    
    for label in labels:
        if label not in label_to_idx:
            continue
        idx = label_to_idx[label]
        
        if "optimized_coordinates" in df.columns:
            coords = df.loc[idx, "optimized_coordinates"]
            if coords is not None and hasattr(coords, '__len__') and len(coords) > 0:
                mol_coords[label] = coords
                max_atoms = max(max_atoms, len(coords))
    
    if not mol_coords:
        return None
    
    rows = []
    for i in range(max_atoms):
        row_data = {"No": i + 1}
        
        for label in labels:
            if label not in mol_coords:
                row_data[f"{label}|Element"] = ""
                row_data[f"{label}|X"] = ""
                row_data[f"{label}|Y"] = ""
                row_data[f"{label}|Z"] = ""
                continue
            
            coords = mol_coords[label]
            if hasattr(coords, 'iloc') and i < len(coords):  # DataFrame
                atom_row = coords.iloc[i]
                row_data[f"{label}|Element"] = atom_row.get("element", "") if hasattr(atom_row, 'get') else ""
                row_data[f"{label}|X"] = f"{atom_row.get('x', 0):.4f}" if hasattr(atom_row, 'get') else ""
                row_data[f"{label}|Y"] = f"{atom_row.get('y', 0):.4f}" if hasattr(atom_row, 'get') else ""
                row_data[f"{label}|Z"] = f"{atom_row.get('z', 0):.4f}" if hasattr(atom_row, 'get') else ""
            elif isinstance(coords, list) and i < len(coords):
                atom = coords[i]
                if isinstance(atom, dict):
                    row_data[f"{label}|Element"] = atom.get("element", "")
                    row_data[f"{label}|X"] = f"{atom.get('x', 0):.4f}"
                    row_data[f"{label}|Y"] = f"{atom.get('y', 0):.4f}"
                    row_data[f"{label}|Z"] = f"{atom.get('z', 0):.4f}"
            else:
                row_data[f"{label}|Element"] = ""
                row_data[f"{label}|X"] = ""
                row_data[f"{label}|Y"] = ""
                row_data[f"{label}|Z"] = ""
        
        rows.append(row_data)
    
    return pd.DataFrame(rows) if rows else None


def build_horizontal_ir_table(df: pd.DataFrame, labels: List[str], label_to_idx: Dict) -> Optional[pd.DataFrame]:
    """Build horizontal IR table: |No| Freq_mol1 | Int_mol1 | Freq_mol2 | Int_mol2 | ... |"""
    
    mol_ir = {}
    max_modes = 0
    
    for label in labels:
        if label not in label_to_idx:
            continue
        idx = label_to_idx[label]
        
        if "ir" in df.columns:
            ir = df.loc[idx, "ir"]
            if ir is not None and hasattr(ir, 'empty') and not ir.empty:
                mol_ir[label] = ir
                max_modes = max(max_modes, len(ir))
    
    if not mol_ir:
        return None
    
    rows = []
    for i in range(min(max_modes, 100)):  # Limit to 100 modes
        row_data = {"No": i + 1}
        
        for label in labels:
            if label not in mol_ir:
                row_data[f"{label}|Freq"] = ""
                row_data[f"{label}|Int"] = ""
                continue
            
            ir = mol_ir[label]
            if i < len(ir):
                freq_col = "freq_cm-1" if "freq_cm-1" in ir.columns else ir.columns[0]
                int_col = "intensity_km/mol" if "intensity_km/mol" in ir.columns else ("intensity" if "intensity" in ir.columns else ir.columns[1] if len(ir.columns) > 1 else None)
                
                row_data[f"{label}|Freq"] = f"{ir.iloc[i][freq_col]:.1f}"
                if int_col:
                    row_data[f"{label}|Int"] = f"{ir.iloc[i][int_col]:.2f}"
                else:
                    row_data[f"{label}|Int"] = ""
            else:
                row_data[f"{label}|Freq"] = ""
                row_data[f"{label}|Int"] = ""
        
        rows.append(row_data)
    
    return pd.DataFrame(rows) if rows else None


def build_horizontal_raman_table(df: pd.DataFrame, labels: List[str], label_to_idx: Dict) -> Optional[pd.DataFrame]:
    """Build horizontal Raman table: |No| Freq_mol1 | Act_mol1 | Freq_mol2 | Act_mol2 | ... |"""
    
    mol_raman = {}
    max_modes = 0
    
    for label in labels:
        if label not in label_to_idx:
            continue
        idx = label_to_idx[label]
        
        if "raman" in df.columns:
            raman = df.loc[idx, "raman"]
            if raman is not None and hasattr(raman, 'empty') and not raman.empty:
                mol_raman[label] = raman
                max_modes = max(max_modes, len(raman))
    
    if not mol_raman:
        return None
    
    rows = []
    for i in range(min(max_modes, 100)):
        row_data = {"No": i + 1}
        
        for label in labels:
            if label not in mol_raman:
                row_data[f"{label}|Freq"] = ""
                row_data[f"{label}|Act"] = ""
                continue
            
            raman = mol_raman[label]
            if i < len(raman):
                freq_col = "freq_cm-1" if "freq_cm-1" in raman.columns else raman.columns[0]
                act_col = "activity" if "activity" in raman.columns else (raman.columns[1] if len(raman.columns) > 1 else None)
                
                row_data[f"{label}|Freq"] = f"{raman.iloc[i][freq_col]:.1f}"
                if act_col:
                    row_data[f"{label}|Act"] = f"{raman.iloc[i][act_col]:.4f}"
                else:
                    row_data[f"{label}|Act"] = ""
            else:
                row_data[f"{label}|Freq"] = ""
                row_data[f"{label}|Act"] = ""
        
        rows.append(row_data)
    
    return pd.DataFrame(rows) if rows else None


def build_horizontal_tddft_table(df: pd.DataFrame, labels: List[str], label_to_idx: Dict) -> Optional[pd.DataFrame]:
    """Build horizontal TDDFT table: |State| eV_mol1 | nm_mol1 | f_mol1 | eV_mol2 | ... |"""
    
    mol_tddft = {}
    max_states = 0
    
    for label in labels:
        if label not in label_to_idx:
            continue
        idx = label_to_idx[label]
        
        if "tddft_states" in df.columns:
            tddft = df.loc[idx, "tddft_states"]
            if tddft is not None and hasattr(tddft, 'empty') and not tddft.empty:
                mol_tddft[label] = tddft
                if "state" in tddft.columns:
                    max_states = max(max_states, tddft["state"].max() if len(tddft) > 0 else 0)
                else:
                    max_states = max(max_states, len(tddft))
    
    if not mol_tddft:
        return None
    
    rows = []
    for i in range(min(int(max_states), 20)):
        state_num = i + 1
        row_data = {"State": f"S{state_num}"}
        
        for label in labels:
            if label not in mol_tddft:
                row_data[f"{label}|eV"] = ""
                row_data[f"{label}|nm"] = ""
                row_data[f"{label}|f"] = ""
                continue
            
            tddft = mol_tddft[label]
            
            # Find the state
            if "state" in tddft.columns:
                state_rows = tddft[tddft["state"] == state_num]
                if len(state_rows) > 0:
                    s = state_rows.iloc[0]
                    energy = s.get("energy_ev", None)
                    if energy is not None:
                        row_data[f"{label}|eV"] = f"{energy:.3f}"
                        row_data[f"{label}|nm"] = f"{1239.84 / energy:.1f}"
                    else:
                        row_data[f"{label}|eV"] = ""
                        row_data[f"{label}|nm"] = ""
                    row_data[f"{label}|f"] = f"{s.get('fosc', 0):.4f}" if "fosc" in tddft.columns else ""
                else:
                    row_data[f"{label}|eV"] = ""
                    row_data[f"{label}|nm"] = ""
                    row_data[f"{label}|f"] = ""
            else:
                row_data[f"{label}|eV"] = ""
                row_data[f"{label}|nm"] = ""
                row_data[f"{label}|f"] = ""
        
        rows.append(row_data)
    
    return pd.DataFrame(rows) if rows else None


def build_horizontal_internal_coords_table(df: pd.DataFrame, labels: List[str], label_to_idx: Dict) -> Optional[pd.DataFrame]:
    """Build horizontal internal coords table from bonds/angles/dihedrals fields.
    
    The DataFrame has separate 'bonds', 'angles', 'dihedrals' columns, each containing a DataFrame
    with columns like: atoms, value (or similar).
    """
    
    # Determine which coordinate types are available
    coord_types = []
    for ct in ["bonds", "angles", "dihedrals"]:
        if ct in df.columns:
            for label in labels:
                if label in label_to_idx:
                    idx = label_to_idx[label]
                    data = df.loc[idx, ct]
                    if data is not None and hasattr(data, 'empty') and not data.empty:
                        coord_types.append(ct)
                        break
    
    if not coord_types:
        return None
    
    # Build combined table with type column
    all_rows = []
    
    for coord_type in ["bonds", "angles", "dihedrals"]:
        if coord_type not in df.columns:
            continue
        
        # Collect all unique definitions across molecules
        all_defs = set()
        mol_data = {}
        
        for label in labels:
            if label not in label_to_idx:
                continue
            idx = label_to_idx[label]
            data = df.loc[idx, coord_type]
            
            if data is None or (hasattr(data, 'empty') and data.empty):
                continue
            
            mol_data[label] = data
            
            # Find definition column - could be 'atoms', 'definition', or first column
            def_col = None
            for col_name in ["atoms", "definition", "coord"]:
                if col_name in data.columns:
                    def_col = col_name
                    break
            if def_col is None and len(data.columns) > 0:
                def_col = data.columns[0]
            
            if def_col:
                all_defs.update(data[def_col].astype(str).tolist())
        
        if not mol_data:
            continue
        
        # Build rows for this coordinate type
        for defn in sorted(all_defs)[:50]:  # Limit per type
            row_data = {"Type": coord_type.capitalize()[0], "Coordinate": defn}  # B/A/D
            
            for label in labels:
                if label not in mol_data:
                    row_data[f"{label}|Value"] = ""
                    continue
                
                data = mol_data[label]
                
                # Find definition column
                def_col = None
                for col_name in ["atoms", "definition", "coord"]:
                    if col_name in data.columns:
                        def_col = col_name
                        break
                if def_col is None and len(data.columns) > 0:
                    def_col = data.columns[0]
                
                # Find value column
                val_col = None
                for col_name in ["value", "Value", "length", "angle", "dihedral"]:
                    if col_name in data.columns:
                        val_col = col_name
                        break
                if val_col is None and len(data.columns) > 1:
                    val_col = data.columns[1]
                
                if def_col and val_col:
                    match = data[data[def_col].astype(str) == defn]
                    if len(match) > 0:
                        row_data[f"{label}|Value"] = f"{match[val_col].iloc[0]:.4f}"
                    else:
                        row_data[f"{label}|Value"] = ""
                else:
                    row_data[f"{label}|Value"] = ""
            
            all_rows.append(row_data)
    
    if not all_rows:
        return None
    
    result_df = pd.DataFrame(all_rows)
    # Add row number
    result_df.insert(0, "No", range(1, len(result_df) + 1))
    
    return result_df


def build_horizontal_energy_table(df: pd.DataFrame, labels: List[str], label_to_idx: Dict) -> Optional[pd.DataFrame]:
    """Build horizontal energy table: |Property| Value_mol1 | Value_mol2 | ... |"""
    
    energy_cols = ["gibbs_Eh", "single_point_Eh", "enthalpy_Eh", "electronic_energy_Eh", 
                   "zpe_Eh", "thermal_correction_Eh", "homo_energy", "lumo_energy"]
    
    rows = []
    for col in energy_cols:
        row_data = {"Property": col.replace("_", " ").title()}
        has_data = False
        
        for label in labels:
            if label not in label_to_idx:
                row_data[f"{label}|Value"] = ""
                continue
            
            idx = label_to_idx[label]
            if col in df.columns:
                val = df.loc[idx, col]
                if val is not None and not pd.isna(val):
                    row_data[f"{label}|Value"] = f"{val:.6f}"
                    has_data = True
                else:
                    row_data[f"{label}|Value"] = ""
            else:
                row_data[f"{label}|Value"] = ""
        
        if has_data:
            rows.append(row_data)
    
    return pd.DataFrame(rows) if rows else None


def render_export_options(table_df: pd.DataFrame, prefix: str):
    """Render export options for a DataFrame."""
    
    with st.expander("📤 Export Options"):
        col1, col2 = st.columns(2)
        
        with col1:
            csv = table_df.to_csv(index=False, encoding='utf-8-sig')
            st.download_button(
                "📥 Download CSV",
                csv,
                f"comparison_{prefix}.csv",
                "text/csv",
                key=f"compare_{prefix}_csv"
            )
        
        with col2:
            try:
                buffer = io.BytesIO()
                with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
                    table_df.to_excel(writer, index=False, sheet_name='Comparison')
                buffer.seek(0)
                st.download_button(
                    "📥 Download Excel",
                    buffer,
                    f"comparison_{prefix}.xlsx",
                    "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    key=f"compare_{prefix}_xlsx"
                )
            except ImportError:
                st.caption("Install openpyxl for Excel export")
