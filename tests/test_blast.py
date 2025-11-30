import unittest
from unittest.mock import patch
import blast_101_search
import blosum as bl
import smith_waterman_p as SW
from collections import defaultdict



class TestExtendDiagonalBlosum62(unittest.TestCase):
    def test_extend_diagonal_blosum(self):
        """
        """
        test_cases = [
            (62, 19, "BLOSUM 62"),
            (45, 24, "BLOSUM 45"),
            (50, 26, "BLOSUM 50"),
            (90, 21, "BLOSUM 90"),
            (80, 20, "BLOSUM 80")
        ]

        for blosum, b_scores, expected in test_cases:
            with patch.object(blast_101_search, 'ddist', new=bl.BLOSUM(blosum)), \
                 patch.object(SW, 'dist', new=bl.BLOSUM(blosum)), \
                 patch.object(blast_101_search, 'word_size', new=3), \
                 patch.dict(blast_101_search.programme_settings.settings["BLAST"], {
                     "min_extension_score": "0",
                     "max_extension_length": "10"
                 }):
                # Sequences
                s0 = "ACD"
                s1 = "ACD"
                # The starting position (0,0)
                pos_s0_s1 = (0, 0)
                # Call extend_diagonal
                score = blast_101_search.extend_diagonal(pos_s0_s1, s0, s1)

                with self.subTest(f"Testing the function for expected score based on {expected}"):
                    with patch('sys.stderr'):
                        self.assertEqual(
                            score,
                            b_scores,
                            f"Expected a perfect diagonal match of {b_scores}, got {score}"
                        )

    def test_choose_left_vs_right_extension(self):
        """
        Test to verify if the function chooses the correct extension (left or right)
        based on the given BLOSUM score.
        """
        test_cases = [
            (62, 13, "BLOSUM 62"),
            (45, 17, "BLOSUM 45"),
            (50, 18, "BLOSUM 50"),
            (90, 14, "BLOSUM 90"),
            (80, 14, "BLOSUM 80")
        ]

        for blosum, b_scores, expected in test_cases:
            with patch.object(blast_101_search, 'ddist', new=bl.BLOSUM(blosum)), \
                 patch.object(SW, 'dist', new=bl.BLOSUM(blosum)), \
                 patch.object(blast_101_search, 'word_size', new=1), \
                 patch.dict(blast_101_search.programme_settings.settings["BLAST"], {
                     "min_extension_score": "0",
                     "max_extension_length": "10"
                 }):

                # Sequences
                s0 = "ACD"
                s1 = "ACC"
                # The starting position (1,1)
                pos_s0_s1 = (1, 1)
                # Call extend_diagonal
                score = blast_101_search.extend_diagonal(pos_s0_s1, s0, s1)

                with self.subTest(f"Testing to see if the function chooses the right extension (left or right) score based on {expected}"):
                    with patch('sys.stderr'):
                        self.assertEqual(
                            score,
                            b_scores,
                            f"Expected the left extension to be chosen"
                        )

    def test_min_extension_score_stop(self):
        """
        Test to see if extension stops when the current score is below the minimum extension score.
        """
        test_cases = [
            (62, 4, "BLOSUM 62"),
            (45, 5, "BLOSUM 45"),
            (50, 5, "BLOSUM 50"),
            (90, 5, "BLOSUM 90"),
            (80, 5, "BLOSUM 80")
        ]

        for blosum, b_scores, expected in test_cases:
            with patch.object(blast_101_search, 'ddist', new=bl.BLOSUM(blosum)), \
                 patch.object(SW, 'dist', new=bl.BLOSUM(blosum)), \
                 patch.object(blast_101_search, 'word_size', new=1), \
                 patch.dict(blast_101_search.programme_settings.settings["BLAST"], {
                     "min_extension_score": "5",
                     "max_extension_length": "10"
                 }):
                # Sequences
                s0 = "ACD"
                s1 = "APD"
                # The starting position (0,0)
                pos_s0_s1 = (0, 0)
                # Call extend_diagonal
                score = blast_101_search.extend_diagonal(pos_s0_s1, s0, s1)

                with self.subTest(f"Testing to see if cscore <= min_extension_score, extension must stop, based on {expected}"):
                    self.assertEqual(
                        score,
                        b_scores,
                        f"Expected the score to match only the initial pairing"
                    )


class TestProcessBlast(unittest.TestCase):

    def test_no_matches_found(self):
        with patch.object(blast_101_search, 'word_size', new=3), \
             patch.object(blast_101_search, 'query_sequence', new=defaultdict(list)), \
             patch.object(blast_101_search, 'qsequence', new="ACGT"), \
             patch.dict(blast_101_search.programme_settings.settings["BLAST"], {
                 "min_extension_score": "0",
                 "max_extension_length": "10",
             }):
            myline_database = "XXXXXX"
            score = blast_101_search.process_blast(myline_database)
            with self.subTest("No matches"):
                self.assertEqual(score, 0, "Expected 0 if no matches found (count_matches < 2).")

    def test_single_match_found(self):
        with patch.object(blast_101_search, 'word_size', new=3), \
             patch.object(blast_101_search, 'query_sequence', new=defaultdict(list, {"GGG": [0]})), \
             patch.object(blast_101_search, 'qsequence', new="GGG"), \
             patch.dict(blast_101_search.programme_settings.settings["BLAST"], {
                 "min_extension_score": "0",
                 "max_extension_length": "10",
             }):
            myline_database = "GGGXXXX"
            score = blast_101_search.process_blast(myline_database)
            with self.subTest("One match found"):
                self.assertEqual(score, 0, "Expected 0 if only one match was found (count_matches < 2).")

    def test_grouping_single_group_returns_zero(self):
        with patch.object(blast_101_search, 'word_size', new=3), \
             patch.object(blast_101_search, 'query_sequence', new=defaultdict(list, {"AAA": [0]})), \
             patch.object(blast_101_search, 'qsequence', new="AAAAA"), \
             patch.dict(blast_101_search.programme_settings.settings["BLAST"], {
                 "min_extension_score": "0",
                 "max_extension_length": "10",
             }):
            myline_database = "AAAAAA"
            score = blast_101_search.process_blast(myline_database)
            with self.subTest("One group formed - single match"):
                self.assertEqual(score, 0, "Expected 0 score because matches form only a single group (res_store length < 2).")

    def test_multiple_matches_no_diagonal(self):
        with patch.object(blast_101_search, 'word_size', new=3), \
             patch.object(blast_101_search, 'query_sequence', new=defaultdict(list, {
                 "ABC": [0],
                 "XYZ": [10]
             })), \
             patch.object(blast_101_search, 'qsequence', new="ABCXXXXXYZ"), \
             patch.dict(blast_101_search.programme_settings.settings["BLAST"], {
                 "min_extension_score": "0",
                 "max_extension_length": "10",
             }):
            myline_database = "ABCXXXXYZ"
            score = blast_101_search.process_blast(myline_database)
            with self.subTest("No diagonal - first check"):
                self.assertEqual(score, 0, "Expected 0 since positions are not in the right order.")

    def test_multiple_diagonal_matches_positive_score(self):
        with patch.object(blast_101_search, 'word_size', new=3), \
             patch.object(blast_101_search, 'query_sequence', new=defaultdict(list, {
                 "ATG": [0],
                 "TGC": [1]
             })), \
             patch.object(blast_101_search, 'qsequence', new="ATGC"), \
             patch.dict(blast_101_search.programme_settings.settings["BLAST"], {
                 "min_extension_score": "0",
                 "max_extension_length": "10",
             }):
            myline_database = "ATGC"
            score = blast_101_search.process_blast(myline_database)
            with self.subTest("Positive diagonal score"):
                self.assertGreater(score, 0, f"Expected a positive bestscore for diagonal match, got {score}")


if __name__ == "__main__":
    unittest.main()
