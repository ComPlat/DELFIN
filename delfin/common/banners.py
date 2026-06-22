# banner_utils.py
# Common banner generation utilities

from typing import List, Optional

from delfin import __version__


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
                          description: str = "Automates ORCA 6.1.1 calculations",
                          version: Optional[str] = None) -> str:
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

        version_text = version if version is not None else f"Version {__version__}"
        content_parts = ["-***-", author, author_name, institution, description, version_text, "-***-"]
        lines = [border]
        lines.extend(center_line(part) for part in content_parts if part)
        lines.append(border)

        return "\n".join(lines)

    @staticmethod
    def create_compact_banner(author: str = "ComPlat",
                            author_name: str = "M. Hartmann",
                            institution: str = "Karlsruhe Institute of Technology (KIT)",
                            description: str = "Automates ORCA 6.1.1 calculations",
                            version: Optional[str] = None,
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

        version_text = version if version is not None else f"Version {__version__}"
        content_parts = ["-***-", author, author_name, institution, description, version_text, "-***-"]
        lines = [border]
        lines.extend(center_line(part) for part in content_parts if part)
        lines.append(border)

        return "\n".join(lines)


def _indent_lines(lines: List[str], spaces: int) -> str:
    if not spaces:
        return "\n".join(lines)
    padding = " " * spaces
    return "\n".join(f"{padding}{line}" for line in lines)



def build_banner(title: str,
                 *,
                 description: str,
                 version: Optional[str] = None,
                 header_indent: int = 0,
                 info_indent: int = 0,
                 author: str = "ComPlat",
                 author_name: str = "M. Hartmann",
                 institution: str = "Karlsruhe Institute of Technology (KIT)") -> str:
    """Build a banner with a title and info block."""

    version_text = version if version is not None else f"Version {__version__}"
    header_lines = BannerGenerator.create_header_banner(title).splitlines()
    info_lines = BannerGenerator.create_info_banner(
        author=author,
        author_name=author_name,
        institution=institution,
        description=description,
        version=version_text,
    ).splitlines()

    header = _indent_lines(header_lines, header_indent)
    info = _indent_lines(info_lines, info_indent)
    return f"{header}\n\n{info}"



def build_standard_banner(*,
                          description: str = "Automates ORCA 6.1.1, xTB 6.7.1 and CREST 3.0.2 runs",
                          header_indent: int = 0,
                          info_indent: int = 0) -> str:
    """Convenience wrapper for the default DELFIN banner."""

    return build_banner(
        "DELFIN",
        description=description,
        header_indent=header_indent,
        info_indent=info_indent,
    )



def build_occupier_banner(*,
                          description: str = "Automates ORCA 6.1.1 calculations",
                          header_indent: int = 0,
                          info_indent: int = 0) -> str:
    """Convenience wrapper for OCCUPIER banner output."""

    return build_banner(
        "OCCUPIER",
        description=description,
        header_indent=header_indent,
        info_indent=info_indent,
    )
