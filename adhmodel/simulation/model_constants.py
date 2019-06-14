# -*- coding: utf-8 -*-
import param
import panel as pn


class ModelConstants(param.Parameterized):
    kinematic_viscosity = param.Number(
        default=0.00001,
        bounds=(0, None),
        softbounds=(0, 0.0001),
        doc='MP MU: Kinematic molecular viscosity. REQUIRED.',
        precedence=1,
    )
    gravity = param.Number(
        default=9.81,
        bounds=(9.8, 32.20),
        doc='MP G: Gravitational acceleration. REQUIRED.',
        precedence=2,
    )
    density = param.Number(
        default=1000.0,
        bounds=(0, None),
        softbounds=(0, 1100.0),
        doc='MP RHO: Density. REQUIRED.',
        precedence=3,
    )
    enable_wetting_drying = param.Boolean(
        default=True,
        doc='MP DTL: Wetting/drying active.',
        precedence=4,
    )
    wet_dry_limit = param.Number(
        default=0,
        bounds=(0, None),
        softbounds=(0, 1.5),
        doc='MP DTL: Wetting and drying limit.',
        precedence=5,
    )
    mannings_unit_constant = param.Number(
        default=1.0,
        bounds=(1.0, 1.486),
        doc='MP MUC: Mannings unit constant.',
        precedence=6,
    )

    def __init__(self):
        super(ModelConstants, self).__init__()
        self._update_enable_wetting_drying()

    @param.depends('enable_wetting_drying', watch=True)
    def _update_enable_wetting_drying(self):
        if self.enable_wetting_drying:
            self.param.wet_dry_limit.precedence = 5
        else:
            self.param.wet_dry_limit.precedence = -1

    def set_not_required(self, value=False):
        self.enable_wetting_drying = value

    def panel(self):
        return pn.Pane(self.param, show_name=False)
