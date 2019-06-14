# -*- coding: utf-8 -*-
import param
import panel as pn
from .material_wind_properties import MaterialWindProperties


class MaterialProperties(param.Parameterized):
    eddy_viscosity_method = param.ObjectSelector(
        default='Constant (EVS)',
        objects=['Constant (EVS)', 'Estimated (EEV)'],
        doc='MP EVS or MP EEV: Method of defining eddy viscosity.',
        precedence=1,
    )
    vxx_eddy_viscosity = param.Number(
        default=0.0,
        bounds=(0, None),
        doc='MP EVS: Vxx term of the kinematic eddy viscosity tensor.',
        precedence=2,
    )
    vxy_eddy_viscosity = param.Number(
        default=0.0,
        bounds=(0, None),
        doc='MP EVS: Vxy term of the kinematic eddy viscosity tensor.',
        precedence=3,
    )
    vyy_eddy_viscosity = param.Number(
        default=0.0,
        bounds=(0, None),
        doc='MP EVS: Vyy term of the kinematic eddy viscosity tensor.',
        precedence=4,
    )
    estimated_eddy_viscosity_weighting_factor = param.Number(
        default=0.5,
        bounds=(0, None),
        softbounds=(0.1, 1.0),
        doc='MP EEV: Weighting factor, method dependent.',
        precedence=5,
    )
    estimated_eddy_viscosity_method = param.ObjectSelector(
        default=1,
        objects={'1 - Isotropic': 1, '2 - Anisotropic': 2, '3 - Smagorinsky': 3, '4 - Stansby': 4},
        doc='MP EEV: Method of defining estimated eddy viscosity.',
        precedence=6,
    )
    coriolis = param.Boolean(
        default=False,
        doc='MP COR: Coriolis parameter active.',
        precedence=7,
    )
    coriolis_latitude = param.Number(
        default=0.0,
        bounds=(-90,90),
        doc='MP COR: Latitude decimal degrees.',
        precedence=8,
    )
    refinement_tolerance = param.Number(
        default=1.0,
        bounds=(0, None),
        doc='MP SRT: Mesh refinement tolerance. REQUIRED.',
        precedence=9,
    )
    max_refinement_level = param.Integer(
        default=0,
        bounds=(0, None),
        doc='MP ML: Maximum mesh refinement level. REQUIRED.',
        precedence=10,
    )
    # TODO TVS, DPL, DPT, RD, TOR - need documentation

    def __init__(self):
        super(MaterialProperties, self).__init__()
        self.transport_properties = {}
        self.wind_properties = MaterialWindProperties()

        self._update_eddy_viscosity_method()
        self._update_coriolis()

    @param.depends('eddy_viscosity_method', watch=True)
    def _update_eddy_viscosity_method(self):
        self.param.vxx_eddy_viscosity.precedence = -1
        self.param.vxy_eddy_viscosity.precedence = -1
        self.param.vyy_eddy_viscosity.precedence = -1
        self.param.estimated_eddy_viscosity_weighting_factor.precedence = -1
        self.param.estimated_eddy_viscosity_method.precedence = -1
        objects = list(self.param.eddy_viscosity_method.get_range())
        if self.eddy_viscosity_method == objects[0]:
            self.param.vxx_eddy_viscosity.precedence = 2
            self.param.vxy_eddy_viscosity.precedence = 3
            self.param.vyy_eddy_viscosity.precedence = 4
        else:
            self.param.estimated_eddy_viscosity_weighting_factor.precedence = 5
            self.param.estimated_eddy_viscosity_method.precedence = 6

    @param.depends('coriolis', watch=True)
    def _update_coriolis(self):
        self.param.coriolis_latitude.precedence = -1
        if self.coriolis:
            self.param.coriolis_latitude.precedence = 8

    def set_not_required(self, value=False):
        self.coriolis = value

    def panel(self):
        return pn.Pane(self.param, show_name=False)


