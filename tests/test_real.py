"""
Full Integration Test - Downloads from HuggingFace

Run: python tests/test_real.py
"""

import sys
from pathlib import Path
import logging

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Setup logging
from src.logger import setup_logging, get_logger
setup_logging(level=logging.INFO)
logger = get_logger("test_real")


def download_sample_data():
    """Download sample data from HuggingFace."""
    from huggingface_hub import snapshot_download
    
    data_dir = Path("./test_data_hf")
    
    if not data_dir.exists():
        logger.info("Downloading sample data from HuggingFace...")
        snapshot_download(
            repo_id="JauharMz/Orca",
            repo_type="dataset",
            local_dir=str(data_dir),
            local_dir_use_symlinks=False
        )
        logger.info(f"Downloaded to: {data_dir}")
    
    return data_dir


def get_sample_files(data_dir, limit=5):
    """Get list of .out files."""
    files = list(data_dir.rglob("*.out"))[:limit]
    logger.info(f"Found {len(files)} sample files")
    return files


def run_all_tests():
    """Run all tests with HuggingFace data."""
    # Download data first
    data_dir = download_sample_data()
    sample_files = get_sample_files(data_dir, limit=10)
    
    if not sample_files:
        logger.error("No sample files found! Check HuggingFace download.")
        return
    
    from src.parser.factory import ParserFactory
    from src.parser.batch import BatchParser
    from src.analysis.hierarchy_detector import HierarchyDetector
    from src.analysis.partition_detector import PartitionDetector
    from src.viz.energy_diagram import EnergyDiagramVisualizer
    from src.export.data_exporter import DataExporter
    
    logger.info("=" * 60)
    logger.info("ORCA VISUALIZATION - FULL INTEGRATION TEST")
    logger.info("=" * 60)
    
    # 1. Parse single file
    logger.info("[1/6] Parsing single file...")
    factory = ParserFactory()
    result = factory.parse(str(sample_files[0]))
    logger.info(f"  File: {sample_files[0].name}")
    logger.info(f"  Atoms: {len(result.geometry.cart_coords) if result.geometry.cart_coords is not None else 0}")
    logger.info(f"  Energy: {result.energy.gibbs_Eh or result.energy.single_point_Eh}")
    
    # 2. Batch parse
    logger.info("[2/6] Batch parsing...")
    parser = BatchParser()
    df = parser.parse_files([str(f) for f in sample_files])
    logger.info(f"  Parsed: {len(df)} files")
    
    # 3. Hierarchy
    logger.info("[3/6] Detecting hierarchy...")
    hierarchy = HierarchyDetector(df).detect()
    logger.info(f"  Roots: {hierarchy.roots[:5] if hierarchy.roots else 'None'}...")
    
    # 4. Partitions
    logger.info("[4/6] Detecting partitions...")
    partitions = PartitionDetector(df).detect()
    logger.info(f"  States: {list(partitions.by_state.keys())}")
    
    # 5. Visualization
    logger.info("[5/6] Creating visualization...")
    fig = EnergyDiagramVisualizer(df).create_figure()
    fig.write_html("test_output.html")
    logger.info(f"  Saved: test_output.html")
    
    # 6. Export
    logger.info("[6/6] Exporting data...")
    exporter = DataExporter(df)
    exporter.to_json("test_export.json")
    exporter.to_csv("test_export.csv")
    logger.info(f"  Saved: test_export.json, test_export.csv")
    
    logger.info("=" * 60)
    logger.info("✅ ALL TESTS PASSED!")
    logger.info("=" * 60)


if __name__ == "__main__":
    run_all_tests()
