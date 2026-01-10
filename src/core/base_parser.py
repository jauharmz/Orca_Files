"""Base parser abstract class."""

from abc import ABC, abstractmethod
from typing import Any, Dict, Optional
import re

from ..logger import get_logger


class BaseParser(ABC):
    """
    Abstract base class for all ORCA parsers.
    
    Each parser is responsible for extracting a specific type of data
    from the ORCA output text.
    """
    
    def __init__(self, text: str, filename: Optional[str] = None):
        """
        Initialize parser.
        
        Args:
            text: Raw ORCA output text
            filename: Optional filename for logging
        """
        self.text = text
        self.filename = filename
        self.logger = get_logger(self.__class__.__name__)
    
    @abstractmethod
    def parse(self) -> Any:
        """
        Parse the text and return extracted data.
        
        Returns:
            Parsed data (type depends on specific parser)
        """
        pass
    
    def _find_pattern(self, pattern: str, flags: int = 0) -> Optional[re.Match]:
        """Find first match of pattern in text."""
        return re.search(pattern, self.text, flags)
    
    def _find_all_patterns(self, pattern: str, flags: int = 0) -> list:
        """Find all matches of pattern in text."""
        return re.findall(pattern, self.text, flags)
    
    def _log_found(self, item: str, count: Optional[int] = None):
        """Log that item was found."""
        if count is not None:
            self.logger.debug(f"Found {count} {item}")
        else:
            self.logger.debug(f"Found {item}")
    
    def _log_not_found(self, item: str):
        """Log that item was not found."""
        self.logger.debug(f"Not found: {item}")
    
    def _parse_float(self, value: str) -> Optional[float]:
        """Safely parse float, handling scientific notation."""
        try:
            return float(value)
        except (ValueError, TypeError):
            return None
