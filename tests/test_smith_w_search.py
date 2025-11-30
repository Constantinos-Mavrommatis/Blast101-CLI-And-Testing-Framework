import os
import sys
import io
import unittest
from unittest.mock import patch
import main

class TestSWIntegration(unittest.TestCase):
    def test_default_values(self):
        """
        """

        test_argv = [
            "progname",
            "--db", "uniprot_bit2.fasta",
            "--query", "logs/SW_test.fasta",
            "--method", "sw"
        ]

        csv_path = os.path.join("logs", "SWsearch.csv")

        if os.path.exists(csv_path):
            os.remove(csv_path)

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
                self.assertEqual(exit_code, 0, "SW pipeline did not exit with code 0.")

            # Check that SWsearch.csv was created and has content.
            with self.subTest("CSV file existence check"):
                self.assertTrue(
                    os.path.exists(csv_path),
                    "Expected logs/SWsearch.csv to be created, but it does not exist."
                )

            if os.path.exists(csv_path):
                with open(csv_path, "r") as f:
                    lines = f.readlines()

                with self.subTest("CSV file content check"):
                    self.assertGreater(
                        len(lines), 3,
                        "Expected logs/SWsearch.csv to have some content, but it's empty or missing lines."
                    )

            if os.path.exists(csv_path):
                os.remove(csv_path)


if __name__ == "__main__":
    unittest.main()
