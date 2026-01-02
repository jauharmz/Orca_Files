"""
Full Integration Test - Downloads from HuggingFace

Run: python tests/test_real.py
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))


def download_sample_data():
    """Download sample data from HuggingFace."""
    from huggingface_hub import snapshot_download
    
    data_dir = Path("./test_data_hf")
    
    if not data_dir.exists():
        print("Downloading sample data from HuggingFace...")
        snapshot_download(
            repo_id="JauharMz/Orca",
            repo_type="dataset",
            local_dir=str(data_dir),
            local_dir_use_symlinks=False
        )
        print(f"Downloaded to: {data_dir}")
    
    return data_dir


def get_sample_files(data_dir, limit=5):
    """Get list of .out files."""
    files = list(data_dir.rglob("*.out"))[:limit]
    print(f"Found {len(files)} sample files")
    return files


def run_all_tests():
    """Run all tests with HuggingFace data."""
    # Download data first
    data_dir = download_sample_data()
    sample_files = get_sample_files(data_dir, limit=10)
    
    if not sample_files:
        print("No sample files found! Check HuggingFace download.")
        return
    
    from src.parser.factory import ParserFactory
    from src.parser.batch import BatchParser
    from src.analysis.hierarchy_detector import HierarchyDetector
    from src.analysis.partition_detector import PartitionDetector
    from src.viz.energy_diagram import EnergyDiagramVisualizer
    from src.export.data_exporter import DataExporter
    
    print("=" * 60)
    print("ORCA VISUALIZATION - FULL INTEGRATION TEST")
    print("=" * 60)
    
    # 1. Parse single file
    print("\n[1/6] Parsing single file...")
    factory = ParserFactory()
    result = factory.parse(str(sample_files[0]))
    print(f"  File: {sample_files[0].name}")
    print(f"  Atoms: {len(result.geometry.cart_coords) if result.geometry.cart_coords is not None else 0}")
    print(f"  Energy: {result.energy.gibbs_Eh or result.energy.single_point_Eh}")
    
    # 2. Batch parse
    print("\n[2/6] Batch parsing...")
    parser = BatchParser()
    df = parser.parse_files([str(f) for f in sample_files])
    print(f"  Parsed: {len(df)} files")
    
    # 3. Hierarchy
    print("\n[3/6] Detecting hierarchy...")
    hierarchy = HierarchyDetector(df).detect()
    print(f"  Roots: {hierarchy.roots[:5] if hierarchy.roots else 'None'}...")
    
    # 4. Partitions
    print("\n[4/6] Detecting partitions...")
    partitions = PartitionDetector(df).detect()
    print(f"  States: {list(partitions.by_state.keys())}")
    
    # 5. Visualization
    print("\n[5/6] Creating visualization...")
    fig = EnergyDiagramVisualizer(df).create_figure()
    fig.write_html("test_output.html")
    print(f"  Saved: test_output.html")
    
    # 6. Export
    print("\n[6/6] Exporting data...")
    exporter = DataExporter(df)
    exporter.to_json("test_export.json")
    exporter.to_csv("test_export.csv")
    print(f"  Saved: test_export.json, test_export.csv")
    
    print("\n" + "=" * 60)
    print("✅ ALL TESTS PASSED!")
    print("=" * 60)


if __name__ == "__main__":
    run_all_tests()
