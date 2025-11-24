"""
BibTeX Citation File Parser (.bibtex)

Parses ORCA BibTeX citation files.
"""

import re
from dataclasses import dataclass, field


@dataclass
class Citation:
    """Single BibTeX citation."""
    key: str = ""
    entry_type: str = ""
    author: str = ""
    title: str = ""
    journal: str = ""
    year: str = ""
    volume: str = ""
    pages: str = ""
    doi: str = ""


@dataclass
class BibtexData:
    """Parsed BibTeX file data."""
    citations: list[Citation] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            'num_citations': len(self.citations),
            'citations': [
                {
                    'key': c.key,
                    'type': c.entry_type,
                    'author': c.author,
                    'title': c.title,
                    'journal': c.journal,
                    'year': c.year,
                    'volume': c.volume,
                    'pages': c.pages,
                    'doi': c.doi
                }
                for c in self.citations
            ]
        }


def parse_bibtex_file(filepath: str) -> BibtexData:
    """Parse BibTeX file."""
    with open(filepath, 'r') as f:
        content = f.read()
    return parse_bibtex_content(content)


def parse_bibtex_content(content: str) -> BibtexData:
    """Parse BibTeX content from string."""
    data = BibtexData()

    # Find all entries
    entries = re.findall(r'@(\w+)\{([^,]+),([^@]+)\}', content, re.DOTALL)

    for entry_type, key, fields in entries:
        citation = Citation(key=key.strip(), entry_type=entry_type.lower())

        # Parse fields
        field_pattern = r'(\w+)\s*=\s*\{([^}]*)\}'
        for field_name, field_value in re.findall(field_pattern, fields):
            field_name = field_name.lower()
            field_value = field_value.strip()

            if field_name == 'author':
                citation.author = field_value
            elif field_name == 'title':
                citation.title = field_value
            elif field_name == 'journal':
                citation.journal = field_value
            elif field_name == 'year':
                citation.year = field_value
            elif field_name == 'volume':
                citation.volume = field_value
            elif field_name == 'pages':
                citation.pages = field_value
            elif field_name == 'doi':
                citation.doi = field_value

        data.citations.append(citation)

    return data


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        data = parse_bibtex_file(sys.argv[1])
        print(f"Citations: {len(data.citations)}")
        for c in data.citations[:3]:
            print(f"  - {c.author.split(',')[0]} et al., {c.journal} ({c.year})")
