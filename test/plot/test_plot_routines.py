# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

"""
Tests of the plot routines.
"""
import unittest

from lindbladmpo.plot_routines import find_index_nearest_time_within_tolerance


class TestFindIndexNearestTimeWithinTolerance(unittest.TestCase):
    """The class testing the plot routines"""

    def test_nearest_time_within_tolerance(self):
        """
        Test find_index_nearest_time_within_tolerance when the nearest
        time entry is within the specified tolerance.
        """
        time_values = [0.1, 0.2, 0.3, 0.4]
        target = 0.25
        tolerance = 0.05
        index = find_index_nearest_time_within_tolerance(time_values, target, tolerance)
        self.assertEqual(index, 1)  # Match at index 1 is the closest

    def test_nearest_time_outside_tolerance(self):
        """
        Test find_index_nearest_time_within_tolerance when the nearest
        time entry is outside the specified tolerance.
        """
        time_values = [0.1, 0.2, 0.3, 0.4]
        target = 0.25
        tolerance = 0.01
        index = find_index_nearest_time_within_tolerance(time_values, target, tolerance)
        self.assertEqual(index, 1)  # Match at index 1 is the closest

    def test_exact_time_entry(self):
        """
        Test find_index_nearest_time_within_tolerance when the target
        time is an exact match to one of the time entries.
        """
        time_values = [0.1, 0.2, 0.3, 0.4]
        target = 0.2
        tolerance = 0.05
        index = find_index_nearest_time_within_tolerance(time_values, target, tolerance)
        self.assertEqual(index, 1)  # Match at index 1


if __name__ == "__main__":
    unittest.main()
