import unittest
from adhmodel.simulation.time_control import TimeControl


class TestIo(unittest.TestCase):

    def test_dependency_time_step_option(self):
        tc = TimeControl()
        # test to ensure the string objects haven't changed
        base_list = ['Steady state solution (TC STD)', 'Time step series (SERIES DT)',
                     'Auto Time Step Find (TC ATF)']
        curr_list = list(tc.param.time_step_option.objects)
        self.assertListEqual(base_list, curr_list, 'param.ObjectSelector objects have changed.')
        # test dependencies on TC STD
        tc.time_step_option = 'Steady state solution (TC STD)'
        self.assertLess(tc.param.max_time_step_size_time_series.precedence, 0)
        self.assertGreater(tc.param.steady_state_min_time_step_size.precedence, 0)
        self.assertGreater(tc.param.steady_state_max_time_step_size.precedence, 0)
        self.assertLess(tc.param.auto_time_step_find_min_time_step_size.precedence, 0)
        self.assertLess(tc.param.auto_time_step_find_max_time_step_size_series.precedence, 0)
        # test dependencies on SERIES DT
        tc.time_step_option = 'Time step series (SERIES DT)'
        self.assertGreater(tc.param.max_time_step_size_time_series.precedence, 0)
        self.assertLess(tc.param.steady_state_min_time_step_size.precedence, 0)
        self.assertLess(tc.param.steady_state_max_time_step_size.precedence, 0)
        self.assertLess(tc.param.auto_time_step_find_min_time_step_size.precedence, 0)
        self.assertLess(tc.param.auto_time_step_find_max_time_step_size_series.precedence, 0)
        # test dependencies on TC ATF
        tc.time_step_option = 'Auto Time Step Find (TC ATF)'
        self.assertLess(tc.param.max_time_step_size_time_series.precedence, 0)
        self.assertLess(tc.param.steady_state_min_time_step_size.precedence, 0)
        self.assertLess(tc.param.steady_state_max_time_step_size.precedence, 0)
        self.assertGreater(tc.param.auto_time_step_find_min_time_step_size.precedence, 0)
        self.assertGreater(tc.param.auto_time_step_find_max_time_step_size_series.precedence, 0)
        # test dependecies on TC STD
        tc.time_step_option = 'Steady state solution (TC STD)'
        self.assertLess(tc.param.max_time_step_size_time_series.precedence, 0)
        self.assertGreater(tc.param.steady_state_min_time_step_size.precedence, 0)
        self.assertGreater(tc.param.steady_state_max_time_step_size.precedence, 0)
        self.assertLess(tc.param.auto_time_step_find_min_time_step_size.precedence, 0)
        self.assertLess(tc.param.auto_time_step_find_max_time_step_size_series.precedence, 0)
