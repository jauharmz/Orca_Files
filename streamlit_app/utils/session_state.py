"""Session state management utilities."""

import streamlit as st
import pandas as pd
from typing import Optional, List


def init_session_state():
    """Initialize session state variables."""
    
    if "df" not in st.session_state:
        st.session_state.df = None
    
    if "selected_molecules" not in st.session_state:
        st.session_state.selected_molecules = []
    
    if "selected_states" not in st.session_state:
        st.session_state.selected_states = []
    
    if "reaction_nodes" not in st.session_state:
        st.session_state.reaction_nodes = []
    
    if "reaction_edges" not in st.session_state:
        st.session_state.reaction_edges = []


def get_df() -> Optional[pd.DataFrame]:
    """Get the main DataFrame."""
    return st.session_state.get("df", None)


def set_df(df: pd.DataFrame):
    """Set the main DataFrame."""
    st.session_state.df = df


def get_selected_molecules() -> List[str]:
    """Get selected molecule IDs."""
    return st.session_state.get("selected_molecules", [])


def get_filtered_df() -> Optional[pd.DataFrame]:
    """Get DataFrame filtered by global selections."""
    df = get_df()
    if df is None:
        return None
    
    filtered = df.copy()
    
    # Filter by molecules
    selected_mols = get_selected_molecules()
    if selected_mols:
        filtered = filtered[filtered["molecule_id"].isin(selected_mols)]
    
    # Filter by states
    selected_states = st.session_state.get("selected_states", [])
    if selected_states and "optimized_state" in filtered.columns:
        filtered = filtered[filtered["optimized_state"].isin(selected_states)]
    
    return filtered
