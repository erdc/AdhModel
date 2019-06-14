import unittest
from adhmodel.simulation.model_constants import ModelConstants


class TestIo(unittest.TestCase):
    def test_dependency_enable_wetting_drying(self):
        mc = ModelConstants()
        mc.enable_wetting_drying = False
        self.assertLess(mc.param.wet_dry_limit.precedence, 0)
        mc.enable_wetting_drying = True
        self.assertGreater(mc.param.wet_dry_limit.precedence, 0)
