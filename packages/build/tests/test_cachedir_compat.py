"""Regression: cachedir Opener.set_type shadows its own ext() helper."""

from __future__ import annotations

import unittest

from omnipath_build.cachedir_compat import patch_cachedir_opener


class TestCachedirCompat(unittest.TestCase):
    def test_missing_ext_does_not_raise(self):
        try:
            from cachedir._open import Opener
        except ImportError:
            self.skipTest("cachedir is not installed")

        patch_cachedir_opener()
        opener = Opener.__new__(Opener)
        opener.ext = None
        opener.path = "download_complexes.php"
        opener.set_type()
        self.assertEqual(opener.ext, "php")
        self.assertEqual(opener.type, "plain")

        opener.ext = None
        opener.path = "signor_complexes"
        opener.set_type()
        self.assertEqual(opener.ext, "txt")
        self.assertEqual(opener.type, "plain")


if __name__ == "__main__":
    unittest.main()
