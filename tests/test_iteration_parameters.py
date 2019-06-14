import unittest
from adhmodel.simulation.iteration_parameters import IterationParameters


class TestIo(unittest.TestCase):
    def test_dependency_non_linear_tolerance_option(self):
        ip = IterationParameters()
        ip.non_linear_tolerance_option = 'Specify residual and incremental (IP NTL & IP ITL)'
        self.assertGreater(ip.param.non_linear_residual_tolerance.precedence, 0)
        self.assertGreater(ip.param.non_linear_incremental_tolerance.precedence, 0)
        ip.non_linear_tolerance_option = 'Specify residual (IP NTL)'
        self.assertGreater(ip.param.non_linear_residual_tolerance.precedence, 0)
        self.assertLess(ip.param.non_linear_incremental_tolerance.precedence, 0)
        ip.non_linear_tolerance_option = 'Specify incremental (IP ITL)'
        self.assertLess(ip.param.non_linear_residual_tolerance.precedence, 0)
        self.assertGreater(ip.param.non_linear_incremental_tolerance.precedence, 0)
        ip.non_linear_tolerance_option = 'Specify residual and incremental (IP NTL & IP ITL)'
        self.assertGreater(ip.param.non_linear_residual_tolerance.precedence, 0)
        self.assertGreater(ip.param.non_linear_incremental_tolerance.precedence, 0)
