# tests/test_cli.py
import sys
import os
import unittest
from unittest.mock import patch
import cli
from cli import validate_blosum, validate_gap, validate_max_alignments, validate_blast_maximum_score, validate_word_size, validate_sw_max_score
import tempfile

class TestCLI(unittest.TestCase):
    def test_default_values(self):
        """
        If the user provides no arguments, the defaults from settings.ini should be used.
        """
        test_argv = ["progname"]
        with patch.object(sys, 'argv', test_argv):
            args = cli.parse_args()
            # Check that defaults match what's in your code:
            self.assertEqual(args.db, cli.default_database)
            self.assertEqual(args.query, cli.default_query)
            self.assertEqual(args.method, "blast")
            self.assertEqual(args.blosum, validate_blosum(cli.default_blosum))
            self.assertEqual(args.gap, validate_gap(cli.default_gap))
            self.assertEqual(args.max_alignments, validate_max_alignments(cli.default_blast_max_alignments))
            self.assertEqual(args.max_sw_scores, validate_sw_max_score(cli.default_sw_max_score))
            self.assertEqual(args.word_size, validate_word_size(cli.default_blast_word_size))
            self.assertEqual(args.max_scores, validate_blast_maximum_score(cli.default_blast_max_score))

    def test_method(self):
        """
        Simulate calling with --method.
        """
        test_cases = [
            (["progname", "--method", "sw"], "sw"),
            (["progname", "--method", "blast"], "blast"),
            (["progname", "--method", "stats"], "stats"),
            (["progname", "--method", "tests"], "tests")
        ]

        for argv, expected in test_cases:
            with self.subTest(argv=argv):
                with patch.object(sys, 'argv', argv):
                    args = cli.parse_args()
                    self.assertEqual(args.method, expected)

        with self.subTest(argv=["progname", "--method", "invalid"]):
            with patch.object(sys, 'argv', ["progname", "--method", "invalid"]):
                with patch('sys.stderr'):
                    with self.assertRaises(SystemExit, msg="Expected parse_args() to exit on invalid method."):
                        cli.parse_args()

    def test_database(self):
        test_cases = [
            (">test_sequence\n****\n", "Database file with invalid characters"),
            ("", "Empty database file"),
            (None, "Nonexistent database file")
        ]

        for content, description in test_cases:
            with self.subTest(description=description):
                if content is not None:
                    tmp_file = tempfile.NamedTemporaryFile(mode="w", delete=False)
                    try:
                        tmp_file.write(content)
                        tmp_file.close()
                        tmp_path = tmp_file.name
                    except Exception:
                        tmp_file.close()
                        os.remove(tmp_file.name)
                        raise
                else:
                    tmp_path = "nonexistent_file.fasta"

                test_argv = ["progname", "--db", tmp_path]
                with patch.object(sys, 'argv', test_argv):
                    with patch('sys.stdout'):
                        with self.assertRaises(SystemExit):
                            cli.parse_args()

                if content is not None:
                    os.remove(tmp_path)


    def test_invalid_blosum(self):
        """
        Test that passing an invalid BLOSUM value raises an error.
        """

        test_argv = [
            ["progname", "--blosum", "999"],
            ["progname", "--blosum", "abc"]
        ]
        for argv in test_argv:
            with self.subTest(argv=argv):
                with patch.object(sys, 'argv', argv):
                    # Because parse_args() calls 'exit' or raises an ArgumentTypeError,
                    # we catch it using assertRaises.
                    with patch('sys.stderr'):
                        with self.assertRaises(SystemExit, msg="Expected parse_args() to exit on invalid BLOSUM."):  # or ArgumentTypeError if you prefer
                            cli.parse_args()

    def test_gap_out_of_range(self):
        """
        If user supplies a gap outside [-20, -1], parse_args should fail.
        """
        test_argv = ["progname", "--gap", "0"]
        with patch.object(sys, 'argv', test_argv):
            with patch('sys.stderr'):
                with self.assertRaises(SystemExit, msg="Expected parse_args() to exit on invalid gap score."):
                    cli.parse_args()

    def test_query_validation(self):
        """
        Check that if user provides a nonexistent query file, we exit with an error.
        """
        test_cases = [
            (">test_sequence\n*****\n", "Database file with invalid characters"),
            ("", "Empty database file"),
            (">test_sequence\nATGCATGGCCGTA\n", "Database file with only DNA characters"),
            (None, "Nonexistent database file")
        ]
        for content, description in test_cases:
            with self.subTest(description=description):
                if content is not None:
                # Create a temporary file with the given content.
                    tmp_file = tempfile.NamedTemporaryFile(mode="w", delete=False)
                    try:
                        tmp_file.write(content)
                        tmp_file.close()
                        tmp_path=tmp_file.name
                    except Exception:
                        tmp_file.close()
                        os.remove(tmp_file.name)
                        raise
                else:
                    tmp_path="nofile.fasta"

                test_argv = ["progname", "--query", tmp_path]
                with patch.object(sys, 'argv', test_argv):
                    with patch('sys.stdout'):
                        with self.assertRaises(SystemExit, msg="Expected parse_args() to exit on invalid query file."):
                            cli.parse_args()

                if content is not None:
                    os.remove(tmp_path)

    def test_sw_max_score(self):
        test_argv = ["progname", "--max_sw_scores", "50000"]
        with patch.object(sys, 'argv', test_argv):
            with patch('sys.stderr'):
                with self.assertRaises(SystemExit, msg="Expected parse_args() to exit on invalid SW max score."):
                    cli.parse_args()

    def test_word_size(self):
        test_argv = ["progname", "--word_size", "100"]
        with patch.object(sys, 'argv', test_argv):
            with patch('sys.stderr'):
                with self.assertRaises(SystemExit, msg="Expected parse_args() to exit on invalid Blast word size."):
                    cli.parse_args()

    def test_blast_maximum_score(self):
        test_argv = ["progname", "--max_scores", "50000"]
        with self.subTest(description=test_argv):
            with patch.object(sys, 'argv', test_argv):
                with patch('sys.stderr'):
                    with self.assertRaises(SystemExit, msg="Expected parse_args() to exit on invalid Blast max score."):
                        cli.parse_args()

        test_argv = ["progname", "--max_scores", "50", "--max_alignments", "500"]
        with self.subTest(description=test_argv):
            with patch.object(sys, 'argv', test_argv):
                with patch('sys.stdout'):
                    with self.assertRaises(SystemExit, msg="Expected parse_args() to exit on invalid Blast max score."):
                        cli.parse_args()

    def test_max_alignments(self):
        test_argv = ["progname", "--max_alignments", "2000"]
        with patch.object(sys, 'argv', test_argv):
            with patch('sys.stderr'):
                with self.assertRaises(SystemExit, msg="Expected parse_args() to exit on invalid Blast max alignments."):
                    cli.parse_args()

    def test_all_arguments(self):
        test_argv = [
            "progname",
            "--db", "uniprot_sprot.fasta",
            "--query", "logs/toyfasta.fasta",
            "--method", "blast",
            "--blosum", "62",
            "--gap", "-10",
            "--max_alignments", "50",
            "--max_scores", "200",
            "--word_size", "4",
            "--max_sw_scores", "500"
        ]
        with patch.object(sys, 'argv', test_argv):
            with patch('sys.stderr'):
                args = cli.parse_args()
                self.assertEqual(args.db, "uniprot_sprot.fasta")
                self.assertEqual(args.query, "pwnaaplhnfgedflqpyvqlqqnfsasdlevnleatreshahfstpqalelflnysvtp")
                self.assertEqual(args.method, "blast")
                self.assertEqual(args.blosum, 62)  # assuming validate_blosum converts to int
                self.assertEqual(args.gap, -10)
                self.assertEqual(args.max_alignments, 50)
                self.assertEqual(args.max_scores, 200)
                self.assertEqual(args.word_size, 4)
                self.assertEqual(args.max_sw_scores, 500)

if __name__ == "__main__":
    unittest.main()
