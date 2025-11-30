import os
import sys
import io
import unittest
from unittest.mock import patch
import main

class TestBlastIntegration(unittest.TestCase):
    def test_default_values(self):
        """
        """

        test_argv = [
            "progname",
            "--db", "uniprot_bit2.fasta",
            "--query", "logs/Blast_test.fasta",
            "--method", "blast"
        ]

        csv_path = os.path.join("logs", "BLsearch.csv")

        if os.path.exists(csv_path):
            os.remove(csv_path)

        with patch.object(sys, 'argv', test_argv):

            # Patch both stdout and stderr to suppress *all* printed output.
            # Instead of printing to the console, any print() calls are written to mock_stdout, mock_stderr.
            mock_stdout = io.StringIO()
            mock_stderr = io.StringIO()
            with patch('sys.stdout', new=mock_stdout), patch('sys.stderr', new=mock_stderr):
                with self.assertRaises(SystemExit) as cm:
                    main.main()

            exit_code = cm.exception.code

            # Check that the pipeline signaled success (exit code=0).
            with self.subTest("Return code check"):
                self.assertEqual(exit_code, 0, "Blast pipeline did not exit with code 0.")

            # Check that BLsearch.csv was created and has content.
            with self.subTest("CSV file existence check"):
                self.assertTrue(
                    os.path.exists(csv_path),
                    "Expected logs/BLsearch.csv to be created, but it does not exist."
                )

            if os.path.exists(csv_path):
                with open(csv_path, "r") as f:
                    lines = f.readlines()

                with self.subTest("CSV file content check"):
                    self.assertGreater(
                        len(lines), 3,
                        "Expected logs/BLsearch.csv to have some content, but it's empty or missing lines."
                    )

            if os.path.exists(csv_path):
                os.remove(csv_path)


if __name__ == "__main__":
    unittest.main()
