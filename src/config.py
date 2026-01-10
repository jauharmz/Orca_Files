"""Global configuration."""

from pathlib import Path


# Paths
PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data"
LEGACY_DIR = PROJECT_ROOT / "legacy"
TEST_DATA_DIR = PROJECT_ROOT / "tests" / "test_data"

# HuggingFace
HF_REPO_ID = "JauharMz/Orca"

# Logging
LOG_LEVEL = "INFO"
LOG_FILE = "orca_viz.log"

# Parser
SUPPORTED_EXTENSIONS = [".out"]

# Visualization
DEFAULT_THEME = "plotly_white"
MOLECULE_STYLES = ["stick", "sphere", "cartoon"]

# Export
EXPORT_FORMATS = ["json", "csv", "parquet", "pkl"]
