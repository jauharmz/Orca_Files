"""Custom exceptions for ORCA parser."""


class ORCAError(Exception):
    """Base exception for ORCA parser."""
    pass


class ORCAFileError(ORCAError):
    """File-related errors (not found, invalid format)."""
    pass


class ORCAParseError(ORCAError):
    """Parsing errors (regex failed, unexpected format)."""
    pass


class ORCAValidationError(ORCAError):
    """Validation errors (invalid data, missing fields)."""
    pass


class ORCAExportError(ORCAError):
    """Export errors (file write failed, invalid format)."""
    pass
