class PepLinkError(Exception):
    """Base exception for PepLink."""


class ValidationError(PepLinkError):
    """Raised when a user input is malformed or inconsistent."""


class UnsupportedPeptideError(PepLinkError):
    """Raised when a peptide is outside PepLink v1's supported scope."""


class MissingDependencyError(PepLinkError):
    """Raised when an optional dependency is required for a feature."""
