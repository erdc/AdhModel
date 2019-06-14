import unittest
from adhmodel.simulation.material_properties import MaterialProperties
from adhmodel.simulation.simulation import BoundaryConditions


class TestIo(unittest.TestCase):
    def test_dependency_eddy_viscosity_method(self):
        mp = MaterialProperties()
        mp.eddy_viscosity_method = 'Constant (EVS)'
        self.assertGreater(mp.param.vxx_eddy_viscosity.precedence, 0)
        self.assertGreater(mp.param.vxy_eddy_viscosity.precedence, 0)
        self.assertGreater(mp.param.vyy_eddy_viscosity.precedence, 0)
        self.assertLess(mp.param.estimated_eddy_viscosity_weighting_factor.precedence, 0)
        self.assertLess(mp.param.estimated_eddy_viscosity_method.precedence, 0)
        mp.eddy_viscosity_method = 'Estimated (EEV)'
        self.assertLess(mp.param.vxx_eddy_viscosity.precedence, 0)
        self.assertLess(mp.param.vxy_eddy_viscosity.precedence, 0)
        self.assertLess(mp.param.vyy_eddy_viscosity.precedence, 0)
        self.assertGreater(mp.param.estimated_eddy_viscosity_weighting_factor.precedence, 0)
        self.assertGreater(mp.param.estimated_eddy_viscosity_method.precedence, 0)
        mp.eddy_viscosity_method = 'Constant (EVS)'
        self.assertGreater(mp.param.vxx_eddy_viscosity.precedence, 0)
        self.assertGreater(mp.param.vxy_eddy_viscosity.precedence, 0)
        self.assertGreater(mp.param.vyy_eddy_viscosity.precedence, 0)
        self.assertLess(mp.param.estimated_eddy_viscosity_weighting_factor.precedence, 0)
        self.assertLess(mp.param.estimated_eddy_viscosity_method.precedence, 0)

    def test_dependency_coriolis(self):
        mp = MaterialProperties()
        mp.coriolis = False
        self.assertFalse(mp.coriolis)
        self.assertLess(mp.param.coriolis_latitude.precedence, 0)
        mp.coriolis = True
        self.assertGreater(mp.param.coriolis_latitude.precedence, 0)

    def test_add_material(self):
        # create default material
        mp = MaterialProperties()
        # create default simulation
        bc_class = BoundaryConditions()
        # add new material with default numbering
        bc_class.add_material(mp)
        self.assertEqual(bc_class.material_properties[2], mp)
        # add new material with specified numbering
        bc_class.add_material(mp, 44)
        self.assertEqual(bc_class.material_properties[44], mp)
