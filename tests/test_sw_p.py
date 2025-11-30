import unittest
from unittest.mock import patch
import smith_waterman_p
import blosum as bl


class TestSmithWaterman(unittest.TestCase):

    def test_single_char_sequences(self):
        """
        """

        test_cases = [
            (62, 4, "BLOSUM 62"),
            (45, 5, "BLOSUM 45"),
            (50, 5, "BLOSUM 50"),
            (90, 5, "BLOSUM 90"),
            (80, 5, "BLOSUM 80")
        ]

        for blosum, b_scores, expected in test_cases:
            smith_waterman_p.dist = bl.BLOSUM(blosum)

            seq1 = "A"
            seq2 = "A"

            with patch.object(smith_waterman_p.args, 'gap', -2):
                score = smith_waterman_p.perform_smith_waterman(
                    seq1, seq2, print_m=False, print_a=False
                )
                with self.subTest(f"Testing the function for expected score based on {expected}"):
                    self.assertEqual(score, b_scores, f"Expected a score of {b_scores} for 'A' vs 'A', got {score}")

    def test_deletion_sequence(self):
        """
        """

        test_cases = [
            (62, 12, "BLOSUM 62"),
            (45, 16, "BLOSUM 45"),
            (50, 17, "BLOSUM 50"),
            (90, 12, "BLOSUM 90"),
            (80, 12, "BLOSUM 80")
        ]

        for blosum, b_scores, expected in test_cases:
            smith_waterman_p.dist = bl.BLOSUM(blosum)

            seq1 = "CPA"
            seq2 = "ACTA"


            with patch.object(smith_waterman_p.args, 'gap', -2):
                score = smith_waterman_p.perform_smith_waterman(
                    seq1, seq2, print_m=False, print_a=False
                )
                with self.subTest(f"Testing the function for a first character deletion scenario on {expected}"):
                    self.assertEqual(score, b_scores, f"Expected a score of {b_scores}, got {score}")

    def test_insertion_sequence(self):
        """
        """

        test_cases = [
            (62, 6, "BLOSUM 62"),
            (45, 8, "BLOSUM 45"),
            (50, 8, "BLOSUM 50"),
            (90, 8, "BLOSUM 90"),
            (80, 8, "BLOSUM 80")
        ]

        for blosum, b_scores, expected in test_cases:
            smith_waterman_p.dist = bl.BLOSUM(blosum)

            seq1 = "AA"
            seq2 = "ACA"

            with patch.object(smith_waterman_p.args, 'gap', -2):
                score = smith_waterman_p.perform_smith_waterman(
                    seq1, seq2, print_m=False, print_a=False
                )
                with self.subTest(f"Testing the function for gap insertion scenario on {expected}"):
                    self.assertEqual(score, b_scores, f"Expected a score of {b_scores}, got {score}")

    def test_multi_char_sequences(self):
        """
        """

        test_cases = [
            (62, 16, "BLOSUM 62"),
            (45, 21, "BLOSUM 45"),
            (50, 22, "BLOSUM 50"),
            (90, 17, "BLOSUM 90"),
            (80, 17, "BLOSUM 80")
        ]

        for blosum, b_scores, expected in test_cases:
            smith_waterman_p.dist = bl.BLOSUM(blosum)

            seq1 = "ACPAA"
            seq2 = "AACTA"

            with patch.object(smith_waterman_p.args, 'gap', -2):
                score = smith_waterman_p.perform_smith_waterman(
                    seq1, seq2, print_m=False, print_a=False
                )
                with self.subTest(f"Testing the function for expected score based on {expected}"):
                    self.assertEqual(score, b_scores, f"Expected a score of {b_scores} for 'ACPAA' vs 'AACTA', got {score}")

    def test_gap_score_logic(self):
        """
        """

        test_cases = [
            (62, 4, "BLOSUM 62"),
            (45, 5, "BLOSUM 45"),
            (50, 5, "BLOSUM 50"),
            (90, 5, "BLOSUM 90"),
            (80, 5, "BLOSUM 80")
        ]

        for blosum, b_scores, expected in test_cases:
            smith_waterman_p.dist = bl.BLOSUM(blosum)

            seq1 = "AA"
            seq2 = "ATA"

            with patch.object(smith_waterman_p.args, 'gap', -5):
                score = smith_waterman_p.perform_smith_waterman(
                    seq1, seq2, print_m=False, print_a=False
                )
                with self.subTest(f"Testing the function for a specific logic based on the gap score on {expected}"):
                    self.assertEqual(score, b_scores, f"Expected a score of {b_scores}, got {score}")

if __name__ == "__main__":
    unittest.main()