import unittest
from adhmodel.simulation.constituent_properties import ConstituentProperties


class TestIo(unittest.TestCase):
    def test_dependency_salinity(self):
        cp = ConstituentProperties()
        cp.salinity = False
        self.assertLess(cp.param.salinity_id.precedence, 0)
        self.assertLess(cp.param.reference_concentration.precedence, 0)
        cp.salinity = True
        self.assertGreater(cp.param.salinity_id.precedence, 0)
        self.assertGreater(cp.param.reference_concentration.precedence, 0)

    def test_dependency_temperature(self):
        cp = ConstituentProperties()
        cp.temperature = False
        self.assertLess(cp.param.temperature_id.precedence, 0)
        self.assertLess(cp.param.reference_temperature.precedence, 0)
        self.assertLess(cp.param.air_water_heat_transfer.precedence, 0)
        self.assertLess(cp.param.short_wave_radiation_series.precedence, 0)
        self.assertLess(cp.param.dew_point_temperature_series.precedence, 0)
        cp.temperature = True
        self.assertGreater(cp.param.temperature_id.precedence, 0)
        self.assertGreater(cp.param.reference_temperature.precedence, 0)
        self.assertGreater(cp.param.air_water_heat_transfer.precedence, 0)
        self.assertGreater(cp.param.short_wave_radiation_series.precedence, 0)
        self.assertGreater(cp.param.dew_point_temperature_series.precedence, 0)

    def test_dependency_vorticity(self):
        cp = ConstituentProperties()
        cp.vorticity = False
        self.assertLess(cp.param.vorticity_id.precedence, 0)
        self.assertLess(cp.param.vorticity_normalization.precedence, 0)
        self.assertLess(cp.param.vorticity_as_term.precedence, 0)
        self.assertLess(cp.param.vorticity_ds_term.precedence, 0)
        cp.vorticity = True
        self.assertGreater(cp.param.vorticity_id.precedence, 0)
        self.assertGreater(cp.param.vorticity_normalization.precedence, 0)
        self.assertGreater(cp.param.vorticity_as_term.precedence, 0)
        self.assertGreater(cp.param.vorticity_ds_term.precedence, 0)

