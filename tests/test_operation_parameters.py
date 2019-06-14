import unittest
from adhmodel.simulation.operation_parameters import OperationParameters


class TestIo(unittest.TestCase):
    def test_dependency_vessel(self):
        op = OperationParameters()
        op.vessel = False
        self.assertLess(op.param.vessel_entrainment.precedence, 0)
        op.vessel = True
        self.assertGreater(op.param.vessel_entrainment.precedence, 0)

    def test_dependency_second_order_temporal_coefficient_active(self):
        op = OperationParameters()
        op.second_order_temporal_coefficient_active = False
        self.assertLess(op.param.second_order_temporal_coefficient.precedence, 0)
        op.second_order_temporal_coefficient_active = True
        self.assertGreater(op.param.second_order_temporal_coefficient.precedence, 0)

    def test_dependency_petrov_galerkin_coefficient_active(self):
        op = OperationParameters()
        op.petrov_galerkin_coefficient_active = False
        self.assertLess(op.param.petrov_galerkin_coefficient.precedence, 0)
        op.petrov_galerkin_coefficient_active = True
        self.assertGreater(op.param.petrov_galerkin_coefficient.precedence, 0)

