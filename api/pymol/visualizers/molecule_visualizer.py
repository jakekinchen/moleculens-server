"""Refactored Molecule Visualizer using Jinja2 templates."""

import logging
import os
from typing import Any, Optional

from jinja2 import Environment, FileSystemLoader

logger = logging.getLogger(__name__)


class MoleculeVisualizer:
    """Molecule visualizer using external templates for better maintainability."""

    def __init__(self, enable_annotations: bool = True):
        """Initialize the molecule visualizer.

        Args:
            enable_annotations: Whether to enable atomic annotations by default
        """
        self.enable_annotations = enable_annotations

        # Set up Jinja2 environment
        template_dir = os.path.join(os.path.dirname(__file__), "..", "templates")
        self.env = Environment(loader=FileSystemLoader(template_dir))

    def _escape_js_string(self, text: str) -> str:
        """Safely escape strings for JavaScript embedding.

        Args:
            text: Text to escape

        Returns:
            Escaped text safe for JavaScript strings
        """
        return (
            text.replace("\\", "\\\\").replace("`", "\\`").replace("$", "\\$").replace("</script>", "<\\/script>")
        )  # Fix HTML injection safety

    def _load_css_styles(self) -> str:
        """Load CSS styles from template file.

        Returns:
            CSS styles as string
        """
        try:
            with open(
                os.path.join(os.path.dirname(__file__), "..", "templates", "molecule_styles.css"),
            ) as f:
                return f.read()
        except Exception as e:
            logger.error(f"Error loading CSS styles: {str(e)}")
            return ""

    def generate_js_code_from_pdb(self, pdb_data: str, name: str = "Molecule") -> str:
        """Generate minimal JavaScript code for embedding in existing Three.js scenes.

        Args:
            pdb_data: PDB format molecular structure data
            name: Name of the molecule

        Returns:
            JavaScript code as string
        """
        try:
            template = self.env.get_template("molecule_viewer.js")

            return template.render(
                pdb_data=self._escape_js_string(pdb_data),
                molecule_name=self._escape_js_string(name),
                molecule_styles=self._escape_js_string(self._load_css_styles()),
                enable_annotations=str(self.enable_annotations).lower(),
            )
        except Exception as e:
            logger.error(f"Error generating JavaScript code: {str(e)}")
            raise ValueError(f"Failed to generate JavaScript code: {str(e)}") from e

    def generate_html_viewer_from_pdb(self, pdb_data: str, name: str = "Molecule") -> str:
        """Generate standalone HTML viewer for a molecule.

        Args:
            pdb_data: PDB format molecular structure data
            name: Name of the molecule

        Returns:
            Complete HTML document as string
        """
        try:
            template = self.env.get_template("viewer.html")

            return template.render(
                pdb_data=self._escape_js_string(pdb_data),
                molecule_name=self._escape_js_string(name),
                molecule_styles=self._load_css_styles(),
                enable_annotations=str(self.enable_annotations).lower(),
            )
        except Exception as e:
            logger.error(f"Error generating HTML viewer: {str(e)}")
            raise ValueError(f"Failed to generate HTML viewer: {str(e)}") from e

    def generate_interactive_html(
        self,
        pdb_data: str,
        title: Optional[str] = None,
        script_data: Optional[dict[str, Any]] = None,
        output_path: Optional[str] = None,
    ) -> str:
        """Generate interactive HTML visualization with optional script data.

        Args:
            pdb_data: PDB format molecular structure data
            title: Title for the molecule visualization
            script_data: Optional script data for interactive animation
            output_path: Optional path to save the HTML file

        Returns:
            HTML content as string if output_path is None, otherwise the path to the generated file
        """
        try:
            # For now, use the simple HTML viewer
            # TODO: Implement interactive features with script_data
            html_content = self.generate_html_viewer_from_pdb(pdb_data, title or "Molecule")

            if output_path:
                with open(output_path, "w") as f:
                    f.write(html_content)
                logger.info(f"Interactive visualization generated successfully: {output_path}")
                return output_path

            return html_content

        except Exception as e:
            logger.error(f"Error generating interactive HTML: {str(e)}")
            raise ValueError(f"Failed to generate interactive HTML: {str(e)}") from e

    def toggle_annotations(self, enable: Optional[bool] = None) -> bool:
        """Toggle or set annotation visibility.

        Args:
            enable: If provided, set annotations to this state. If None, toggle current state.

        Returns:
            New annotation state
        """
        if enable is not None:
            self.enable_annotations = enable
        else:
            self.enable_annotations = not self.enable_annotations

        return self.enable_annotations
