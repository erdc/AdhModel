# -*- coding: utf-8 -*-
import param
import panel as pn


class OperationParameters(param.Parameterized):
    physics = param.ObjectSelector(
        default='SW2',
        objects=['SW2', 'SW3'],
        doc='OP SW2 or OP SW3: Simulation physics. REQUIRED.',
        precedence=1,
    )
    incremental_memory = param.Integer(
        default=40,
        bounds=(1, None),
        softbounds=(20, 200),
        doc='OP INC: Incremental memory block size. REQUIRED.',
        precedence=2,
    )
    blocks_per_processor = param.Integer(
        default=1,
        bounds=(1, None),
        softbounds=(1, 10),
        doc='OP BLK: Number of preconditioning blocks per processor. REQUIRED.',
        precedence=4,
    )
    preconditioner_type = param.ObjectSelector(
        default=1,
        objects={'0 - none': 0, '1 - one level Additive Schwarz': 1,
                 '2 - two level Additive Schwarz': 2, '3 - two level Hybrid': 3},
        doc='OP PRE: Solver pre-conditioner type. REQUIRED.',
        precedence=5,
    )
    transport = param.Integer(
        default=0,
        bounds=(0, None),
        doc='OP TRN: Number of transport constituents.',
        precedence=-1
    )
    vessel = param.Boolean(
        default=True,
        doc='OP BT: Vessel movement active.',
        precedence=6,
    )
    vessel_entrainment = param.Boolean(
        default=False,
        doc='OP BTS: Vessel entrainment active.',
        precedence=-1,
    )
    second_order_temporal_coefficient_active = param.Boolean(
        default=True,
        doc='OP TEM: The second order temporal coefficient will be specified.',
        precedence=8,
    )
    second_order_temporal_coefficient = param.Number(
        default=0,
        bounds=(0, 1),
        doc='OP TEM: Second order temporal coefficient.',
        precedence=8.1,
    )
    petrov_galerkin_coefficient_active = param.Boolean(
        default=True,
        doc="OP TPG: Petrov-Galerkin coefficient active.",
        precedence=9,
    )
    petrov_galerkin_coefficient = param.Number(
        default=0.0,
        bounds=(0.0, 0.5),
        doc="OP TPG: Petrov-Galerkin coefficient.",
        precedence=9.1,
    )
    velocity_gradient = param.Boolean(
        default=False,
        doc="OP NF2: Velocity gradient active.",
        precedence=10,
    )
    wind = param.Boolean(
        default=False,
        doc="OP WND: Wind stress active.",
        precedence=11,
    )
    wave = param.Boolean(
        default=False,
        doc="OP WAV: Short wave stress active.",
        precedence=12,
    )
    dam = param.Boolean(
        default=False,
        doc="OP DAM: Dam break stabilization active.",
        precedence=13,
    )
    diffusive_wave = param.Boolean(
        default=False,
        doc="OP DIF: Diffusive wave solver active.",
        precedence=14,
    )

    def __init__(self):
        super(OperationParameters, self).__init__()
        self._update_vessel_entrainment()
        self._update_second_order_temporal_coefficient_active()
        self._update_petrov_galerkin_coefficient_active()

    # vessel_entrainment requires vessel
    @param.depends('vessel', watch=True)
    def _update_vessel_entrainment(self):
        self.param.vessel_entrainment.precedence = -1
        if self.vessel:
            self.param.vessel_entrainment.precedence = 7

    @param.depends('second_order_temporal_coefficient_active', watch=True)
    def _update_second_order_temporal_coefficient_active(self):
        self.param.second_order_temporal_coefficient.precedence = -1
        if self.second_order_temporal_coefficient_active:
            self.param.second_order_temporal_coefficient.precedence = 8.1

    @param.depends('petrov_galerkin_coefficient_active', watch=True)
    def _update_petrov_galerkin_coefficient_active(self):
        self.param.petrov_galerkin_coefficient.precedence = -1
        if self.petrov_galerkin_coefficient_active:
            self.param.petrov_galerkin_coefficient.precedence = 9.1

    def set_not_required(self, value=False):
        self.vessel = value
        self.vessel_entrainment = value
        self.second_order_temporal_coefficient_active = value
        self.petrov_galerkin_coefficient_active = value
        self.velocity_gradient = value
        self.wind = value
        self.wave = value
        self.dam = value
        self.diffusive_wave = value

    def panel(self):
        return pn.Pane(self.param, show_name=False)
