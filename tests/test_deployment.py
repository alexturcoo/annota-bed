import io
import os
import tempfile
import unittest
from unittest.mock import patch

import api.index as api


class DeploymentTests(unittest.TestCase):
    def test_csv_is_included_in_response(self):
        with api.app.app_context():
            path = api._write_csv([])
            try:
                response = api._response_payload([], [], {}, {}, path, "hg38", "built-in")
                self.assertTrue(response.get_json()["csv_content"].startswith("region_id,"))
            finally:
                os.remove(path)

    def test_prefixed_routes_reach_flask(self):
        client = api.app.test_client()
        self.assertEqual(client.get("/api/health").status_code, 200)
        self.assertEqual(client.get("/api/references").status_code, 200)

    def test_serverless_missing_dependency_returns_error(self):
        with patch.object(api, "IS_VERCEL", True), patch.object(api, "HAS_PYSAM", False):
            response = api.app.test_client().post("/api/annotate")
        self.assertEqual(response.status_code, 503)
        self.assertNotIn("rows", response.get_json())

    def test_assets_download_once_and_skip_lfs_pointers(self):
        with tempfile.TemporaryDirectory() as root:
            source = os.path.join(root, "source")
            cache = os.path.join(root, "cache")
            os.makedirs(source)
            with open(os.path.join(source, "asset.bed"), "wb") as handle:
                handle.write(b"version https://git-lfs.github.com/spec/v1\n")
            data = b"chr1\t10\t20\tID1\tID2\tPLS\n"
            remote = io.BytesIO(data)
            remote.headers = {"Content-Length": str(len(data))}
            with patch.object(api, "BASE_DIR", source), patch.object(api, "IS_VERCEL", True), \
                    patch.object(api, "ASSET_CACHE_DIR", cache), patch.object(api, "urlopen", return_value=remote) as download:
                path = api._ensure_asset("asset.bed")
                self.assertEqual(api._ensure_asset("asset.bed"), path)
                download.assert_called_once()
                with open(path, "rb") as handle:
                    self.assertEqual(handle.read(), data)

    def test_incomplete_download_is_not_cached(self):
        with tempfile.TemporaryDirectory() as root:
            remote = io.BytesIO(b"partial")
            remote.headers = {"Content-Length": "100"}
            with patch.object(api, "BASE_DIR", root), patch.object(api, "IS_VERCEL", True), \
                    patch.object(api, "ASSET_CACHE_DIR", os.path.join(root, "cache")), \
                    patch.object(api, "urlopen", return_value=remote):
                with self.assertRaises(ValueError):
                    api._ensure_asset("asset.bed")
            self.assertFalse(os.path.exists(os.path.join(root, "cache", "asset.bed")))
            self.assertFalse(os.path.exists(os.path.join(root, "cache", "asset.bed.part")))


if __name__ == "__main__":
    unittest.main()
