import os
import sys
import io
import unittest
from unittest.mock import patch
import main
import calc_bit_and_evalues


class TestStatsIntegration(unittest.TestCase):
    def test_default_values(self):
        """
        """

        test_argv = [
            "progname",
            "--method", "stats"
        ]
        with patch.dict(calc_bit_and_evalues.ps.settings["BUILD_EXPECT"], {
                     "show_fit": "False",
                 }):
            # Patch sys.argv so that main.parse_args() sees these arguments.
            with patch.object(sys, 'argv', test_argv):

                mock_stdout = io.StringIO()
                mock_stderr = io.StringIO()
                with patch('sys.stdout', new=mock_stdout), patch('sys.stderr', new=mock_stderr):
                    with self.assertRaises(SystemExit) as cm:
                        main.main()

                exit_code = cm.exception.code

                # Check that the pipeline signaled success (exit code=0).
                with self.subTest("Return code check"):
                    self.assertEqual(exit_code, 0, "Stats pipeline did not exit with code 0.")

if __name__ == "__main__":
    unittest.main()
