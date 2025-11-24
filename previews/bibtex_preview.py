"""BibTeX citation file preview."""

import streamlit as st
import pandas as pd
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from parsers.bibtex_parser import parse_bibtex_file, parse_bibtex_content


def bibtex_preview_page():
    st.title("BibTeX Citation Viewer")

    uploaded_file = st.file_uploader("Upload BibTeX file", type=['bibtex', 'bib'])
    use_sample = st.checkbox("Use sample file")

    if use_sample:
        filepath = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'p1xs0.bibtex')
        data = parse_bibtex_file(filepath)
    elif uploaded_file:
        content = uploaded_file.read().decode('utf-8')
        data = parse_bibtex_content(content)
    else:
        st.info("Please upload a BibTeX file or select a sample file.")
        return

    st.metric("Number of Citations", len(data.citations))

    st.header("Citations")
    for i, c in enumerate(data.citations):
        with st.expander(f"{i+1}. {c.author.split(',')[0] if c.author else 'Unknown'} ({c.year})"):
            st.write(f"**Title:** {c.title}")
            st.write(f"**Authors:** {c.author}")
            st.write(f"**Journal:** {c.journal} {c.volume}, {c.pages}")
            if c.doi:
                st.write(f"**DOI:** [{c.doi}](https://doi.org/{c.doi})")

    if st.button("Export as JSON"):
        import json
        st.download_button("Download JSON", json.dumps(data.to_dict(), indent=2), "citations.json", "application/json")


if __name__ == '__main__':
    bibtex_preview_page()
