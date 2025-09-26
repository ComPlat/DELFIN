# banner_utils.py
# Common banner generation utilities

from typing import Dict, List, Optional


class BannerGenerator:
    """Utility class for generating consistent banners across DELFIN modules."""

    # Standard banner configuration
    BANNER_WIDTH = 61  # Total width of banner including # symbols
    INNER_WIDTH = 59   # Width available for content (excluding # symbols)

    @staticmethod
    def create_header_banner(title: str, subtitle: str = "") -> str:
        """Create a standard header banner.

        Args:
            title: Main title (e.g., "DELFIN", "OCCUPIER")
            subtitle: Optional subtitle

        Returns:
            Formatted banner string
        """
        lines = []

        # Title with asterisks
        title_line = f"*{title.center(len(title) + 4)}*"
        asterisk_line = "*" * len(title_line)

        lines.append(" " * ((BannerGenerator.BANNER_WIDTH - len(asterisk_line)) // 2) + asterisk_line)
        lines.append(" " * ((BannerGenerator.BANNER_WIDTH - len(title_line)) // 2) + title_line)
        lines.append(" " * ((BannerGenerator.BANNER_WIDTH - len(asterisk_line)) // 2) + asterisk_line)

        if subtitle:
            lines.append("")
            lines.append(subtitle)

        return "\n".join(lines)

    @staticmethod
    def create_info_banner(author: str = "ComPlat",
                          author_name: str = "M. Hartmann",
                          institution: str = "Karlsruhe Institute of Technology (KIT)",
                          description: str = "Automates ORCA 6.1.0 calculations",
                          version: str = "Version 1.0.2") -> str:
        """Create a standard info banner with author and institution info.

        Args:
            author: Author/organization name
            institution: Institution name
            description: Brief description
            version: Version string

        Returns:
            Formatted banner string
        """
        border = "#" * BannerGenerator.BANNER_WIDTH

        # Calculate centering for each line
        def center_line(text: str) -> str:
            padding = (BannerGenerator.INNER_WIDTH - len(text)) // 2
            return f"#{' ' * padding}{text}{' ' * (BannerGenerator.INNER_WIDTH - len(text) - padding)}#"

        lines = [
            border,
            center_line("-***-"),
            center_line(author),
            center_line(author_name),
            center_line(institution),
            center_line(description),
            center_line(version),
            center_line("-***-"),
            border
        ]

        return "\n".join(lines)

    @staticmethod
    def create_compact_banner(author: str = "ComPlat",
                            author_name: str = "M. Hartmann",
                            institution: str = "Karlsruhe Institute of Technology (KIT)",
                            description: str = "Automates ORCA 6.1.0 calculations",
                            version: str = "Version 1.0.2",
                            width: int = 49) -> str:
        """Create a compact banner for smaller spaces.

        Args:
            author: Author/organization name
            institution: Institution name
            description: Brief description
            version: Version string
            width: Total banner width

        Returns:
            Formatted banner string
        """
        border = "#" * width
        inner_width = width - 2

        def center_line(text: str) -> str:
            padding = (inner_width - len(text)) // 2
            return f"#{' ' * padding}{text}{' ' * (inner_width - len(text) - padding)}#"

        lines = [
            border,
            center_line("-***-"),
            center_line(author),
            center_line(author_name),
            center_line(institution),
            center_line(description),
            center_line(version),
            center_line("-***-"),
            border
        ]

        return "\n".join(lines)

    @staticmethod
    def create_method_info_block(method_parts: List[str],
                               metals: List[str] = None,
                               metal_basis: str = "",
                               charge: int = 0,
                               multiplicity: int = 1) -> str:
        """Create a formatted method information block.

        Args:
            method_parts: List of method components (functional, basis set, etc.)
            metals: List of transition metals found
            metal_basis: Metal basis set if different
            charge: System charge
            multiplicity: System multiplicity

        Returns:
            Formatted method info string
        """
        lines = []

        # Method line
        method_line = "Method: " + " ".join(part for part in method_parts if part)
        lines.append(method_line)

        # Metal basis if present
        if metals and metal_basis:
            metal_line = f"        {', '.join(metals)} {metal_basis}"
            lines.append(metal_line)

        lines.append("")
        lines.append(f"Charge: {charge}")
        lines.append("-" * 13)

        if multiplicity > 1:
            lines.append(f"Total spin: {(multiplicity - 1) / 2:.1f}")
            lines.append(f"Multiplicity: {multiplicity}")
        else:
            lines.append("Total spin: 0.0")
            lines.append("Multiplicity: 1 (Closed Shell)")

        return "\n".join(lines)

    @staticmethod
    def wrap_in_banner(content: str, title: str = "", width: int = 80) -> str:
        """Wrap content in a simple banner.

        Args:
            content: Content to wrap
            title: Optional title
            width: Banner width

        Returns:
            Content wrapped in banner
        """
        border = "=" * width
        lines = [border]

        if title:
            lines.append(f" {title} ".center(width))
            lines.append(border)

        # Add content with proper indentation
        for line in content.split('\n'):
            if line.strip():
                lines.append(f" {line}")
            else:
                lines.append("")

        lines.append(border)
        return "\n".join(lines)