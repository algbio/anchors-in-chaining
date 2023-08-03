#!/usr/bin/python3
import unittest
from scripts.summary import *


class RunTest(unittest.TestCase):
    def setUp(self) -> None:
        self.chain = [(0, 0, 3), (3, 3, 3), (6, 6, 3)]
        self.chain2 = [(0, 0, 3), (6, 6, 3)]
        self.chain3 = [(9, 9, 3), (11, 11, 3)]
        self.overlap_chain = [(0, 0, 3), (1, 1, 3), (2, 2, 3)]
        self.overlap_chain2 = [(0, 0, 3), (1, 1, 3), (9, 9, 3)]
        self.overlap_chain3 = [(0, 0, 3), (0, 0, 3), (1, 9, 3)]
        self.overlength_chain = [(0, 0, 10)]
        self.test_chain = [(1, 1, 3), (2, 2, 6), (9, 4, 5),
                           (9, 9, 6), (15, 18, 4)]
        self.read_start_pos, self.read_end_pos = 0, 8
        self.real_chain = [(707792, 140, 55), (708062, 415, 28), (708121, 477, 81), (708339, 703, 40), (708386, 749, 62), (708535, 901, 45), (708581, 946, 38), (708621, 987, 58), (708789, 1160, 45), (708843, 1217, 41), (708985, 1365, 104), (709090, 1472, 45), (709135, 1518, 47), (709262, 1654, 30), (709311, 1703, 44), (709357, 1751, 42), (709788, 2207, 40), (709842, 2263, 37), (709879, 2301, 51), (709944, 2368, 49), (710036, 2464, 30), (710138, 2575, 58), (710252, 2694, 39), (710350, 2799, 91), (710479, 2931, 53), (710535, 2988, 40), (710590, 3043, 52), (710643, 3095, 45), (710758, 3214, 35), (710983, 3451, 57), (711208, 3681, 51), (711375, 3853, 36), (711532, 4018, 37), (711611, 4103, 85), (
            711724, 4216, 41), (711764, 4257, 37), (711838, 4334, 31), (711898, 4396, 35), (712037, 4544, 56), (712111, 4622, 38), (712147, 4659, 31), (712178, 4691, 57), (712424, 4943, 35), (712598, 5128, 71), (712688, 5222, 46), (712776, 5313, 44), (712863, 5407, 30), (713143, 5693, 46), (713244, 5800, 46), (713293, 5849, 56), (713346, 5903, 67), (713545, 6113, 49), (713629, 6200, 32), (713723, 6296, 28), (713771, 6346, 31), (713894, 6477, 51), (713944, 6528, 35), (714036, 6625, 28), (714199, 6800, 38), (714252, 6854, 37), (714289, 6892, 125), (714420, 7025, 128), (714663, 7271, 60), (714724, 7331, 47), (714779, 7388, 68), (714868, 7477, 80), (714948, 7558, 37), (714985, 7596, 78)]
        self.real_start_pos = 707658
        self.real_length = 7712

    def test_intersection(self):
        
        read_start_pos2, read_end_pos2 = 10, 18
        intersection_size = intersection(
            self.chain, self.read_start_pos, self.read_end_pos)
        intersection_size2 = intersection(
            self.chain2, self.read_start_pos, self.read_end_pos)
        intersection_size_empty = intersection(
            self.chain, read_start_pos2, read_end_pos2)
        intersection_size_empty2 = intersection(
            self.chain3, self.read_start_pos, self.read_end_pos)
        self.assertEqual(9, intersection_size)
        self.assertEqual(6, intersection_size2)
        self.assertEqual(0, intersection_size_empty)
        self.assertEqual(0, intersection_size_empty2)

        overlap_intersection = intersection(
            self.overlap_chain, self.read_start_pos, self.read_end_pos)
        self.assertEqual(5, overlap_intersection)

        overlap_intersection2 = intersection(
            self.overlap_chain2, self.read_start_pos, self.read_end_pos)
        self.assertEqual(4, overlap_intersection2)

        overlap_intersection3 = intersection(
            self.overlap_chain3, self.read_start_pos, self.read_end_pos)
        self.assertEqual(4, overlap_intersection3)

        overlap_intersection3 = intersection(
            self.overlength_chain, self.read_start_pos, self.read_end_pos)
        self.assertEqual(9, overlap_intersection3)

        overlap_intersection4 = intersection(
            self.chain3, read_start_pos2, read_end_pos2)
        self.assertEqual(5, overlap_intersection4)

        test_start_pos, test_end_pos = 2, 14
        test_intersection = intersection(
            self.test_chain, test_start_pos, test_end_pos)
        self.assertEqual(12, test_intersection)
       
        real_intersection = intersection(self.real_chain, self.real_start_pos, self.real_length)
        print(real_intersection)

    def test_coverage(self):
        coverage1 = coverage(self.chain, 1)
        self.assertEqual(9, coverage1)
        coverage2 = coverage(self.chain2, 1)
        self.assertEqual(6, coverage2)
        coverage3 = coverage(self.chain3, 1)
        self.assertEqual(coverage3, 6)
        coverage4 = coverage(self.overlap_chain, 1)
        self.assertEqual(5, coverage4)
        coverage5 = coverage(self.overlap_chain2, 1)
        self.assertEqual(7, coverage5)
        coverage6 = coverage(self.overlap_chain3, 1)
        self.assertEqual(6, coverage6)
        coverage7 = coverage(self.overlength_chain, 1)
        self.assertEqual(10, coverage7)
        coverage8 = coverage(self.test_chain, 1)
        self.assertEqual(18, coverage8)

    def test_jaccard_index(self):
        jaccard1 = jaccard_index(
            self.read_start_pos, self.read_end_pos, self.chain)
        self.assertEqual(1, jaccard1)
        jaccard2 = jaccard_index(
            self.read_start_pos, self.read_end_pos, self.overlap_chain)
        self.assertEqual(5/9, jaccard2)
        jaccard3 = jaccard_index(
            self.read_start_pos, self.read_end_pos, self.test_chain)
        self.assertEqual(7/(18+9-7), jaccard3)


if __name__ == '__main__':
    unittest.main()
