"""Verify the linked-GEO operations are exposed by the data-normalizer MCP server."""

from __future__ import annotations

import asyncio
import sys
import tempfile
import unittest
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
MCP_ROOT = REPOSITORY_ROOT / "mcp"
if str(MCP_ROOT) not in sys.path:
    sys.path.insert(0, str(MCP_ROOT))

from data_normalizer.server import create_server  # noqa: E402


class DataNormalizerServerTests(unittest.TestCase):
    def test_exposes_linked_geo_tools(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            server = create_server(
                parsed_root=root / "parsed",
                data_root=root / "data",
            )
            tools = asyncio.run(server.list_tools())

        names = {tool.name for tool in tools}
        self.assertIn("fetch_linked_geo", names)
        self.assertIn("parse_linked_geo_matrices", names)
        self.assertIn("fetch_planned_geo", names)
        self.assertIn("parse_geo_matrices", names)


if __name__ == "__main__":
    unittest.main()
