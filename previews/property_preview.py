"""
Property file preview using Streamlit.
"""

import streamlit as st
import plotly.express as px
import pandas as pd
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from parsers.property_parser import parse_property_file, parse_property_content


def property_preview_page():
    """Main Streamlit page for property file preview."""
    st.title("Property File Viewer")

    uploaded_file = st.file_uploader("Upload property file", type=['txt'])
    use_sample = st.checkbox("Use sample file")

    if use_sample:
        filepath = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'p1xs0.property.txt'
        )
        data = parse_property_file(filepath)
    elif uploaded_file:
        content = uploaded_file.read().decode('utf-8')
        data = parse_property_content(content)
    else:
        st.info("Please upload a property.txt file or select a sample file.")
        return

    # Status
    st.header("Calculation Status")
    col1, col2 = st.columns(2)
    with col1:
        st.metric("Status", data.status)
    with col2:
        st.metric("Number of Atoms", data.num_atoms)

    # Charges comparison
    if data.mulliken_charges and data.loewdin_charges:
        st.header("Atomic Charges Comparison")

        charge_data = []
        for i in range(len(data.mulliken_charges)):
            charge_data.append({
                'Atom': i,
                'Mulliken': data.mulliken_charges[i],
                'Loewdin': data.loewdin_charges[i]
            })

        df = pd.DataFrame(charge_data)

        fig = px.bar(
            df.melt(id_vars=['Atom'], var_name='Method', value_name='Charge'),
            x='Atom', y='Charge', color='Method',
            barmode='group',
            title="Mulliken vs Loewdin Charges"
        )
        fig.update_layout(template="plotly_white")
        st.plotly_chart(fig, use_container_width=True)

        # Table
        st.dataframe(df, use_container_width=True)

    # Export
    if st.button("Export as JSON"):
        import json
        st.download_button(
            "Download JSON",
            json.dumps(data.to_dict(), indent=2),
            "property.json",
            "application/json"
        )


if __name__ == '__main__':
    property_preview_page()
