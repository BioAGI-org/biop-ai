# This code is part of the biop-ai distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Offline tests of Bio.TogoWS; the online tests are in test_TogoWS.py."""

import unittest

import support


class TogoWSTimeout(unittest.TestCase):
    def test_timeout_get(self):
        """A silent server must not block a TogoWS GET request forever."""
        support.assert_client_times_out(
            self,
            """
            from Bio import TogoWS
            TogoWS.timeout = 0.5
            TogoWS._open(url)
            """,
        )

    def test_timeout_post(self):
        """A silent server must not block a TogoWS POST request forever."""
        support.assert_client_times_out(
            self,
            """
            from Bio import TogoWS
            TogoWS.timeout = 0.5
            TogoWS._open(url, post="ACGT")
            """,
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
