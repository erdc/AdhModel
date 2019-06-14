# -*- coding: utf-8 -*-
import param
import panel as pn
import pandas as pd


class ConstituentProperties(param.Parameterized):
    salinity = param.Boolean(
        default=False,
        doc='CN SAL: Salinity baroclinic transport active.',
        precedence=1,
    )
    salinity_id = param.Integer(
        default=1,
        bounds=(1, None),
        doc='CN SAL: Constituent ID number.',
        precedence=1.1,
    )
    reference_concentration = param.Number(
        default=35.0,
        bounds=(0, None),
        doc='CN SAL: Reference concentration (ppt).',
        precedence=2,
    )
    temperature = param.Boolean(
        default=False,
        doc='CN TMP: Temperature baroclinic transport active.',
        precedence=3,
    )
    temperature_id = param.Integer(
        default=1,
        bounds=(1, None),
        doc='CN TMP: Constituent ID number.',
        precedence=3.1,
    )
    reference_temperature = param.Number(
        default=12.0,
        softbounds=(-5.0,40),
        doc='CN TMP: Reference temperature degrees C.',
        precedence=4,
    )
    air_water_heat_transfer = param.Boolean(
        default=False,
        doc='CN TMP: Air/water heat transport active.',
        precedence=5,
    )
    short_wave_radiation_series = param.Integer(
        default=0,
        bounds=(0, None),
        doc='DB RAD: Time series id for short wave radiation (w/m^2).',
        precedence=5.1,
    )
    dew_point_temperature_series = param.Integer(
        default=0,
        bounds=(0, None),
        doc='DB RAD: Time series for dew point temperature in degrees C.',
        precedence=5.2,
    )
    vorticity = param.Boolean(
        default=False,
        doc='CN VOR: Vorticity bendway correction active.',
        precedence=6,
    )
    vorticity_id = param.Integer(
        default=1,
        bounds=(1, None),
        doc='CN VOR: Constituent ID number.',
        precedence=6.1,
    )
    vorticity_normalization = param.Number(
        default=1.0,
        bounds=(0, None),
        doc='CN VOR: Normalization factor.',
        precedence=7,
    )
    vorticity_as_term = param.Number(
        default=5.0,
        bounds=(0, None),
        doc='CN VOR: As empirical term.',
        precedence=8,
    )
    vorticity_ds_term = param.Number(
        default=0.5,
        bounds=(0, None),
        doc='CN VOR: Ds empirical term.',
        precedence=9,
    )
    general_constituents = param.DataFrame(
        default=pd.DataFrame({'ID': [], 'CONC': []}),
        doc='CN CON: CARD ID CONCENTRATION - General constituents, '
            'CARD: File card (CN CON), '
            'ID: ID of the constituent, '
            'CONCENTRATION: Concentration of the constituent.',
        precedence=10,
    )
    sand = param.DataFrame(
        default=pd.DataFrame({'ID': [], 'C_0': [], 'C_1': [], 'C_2': [], 'C_3': []}),
        doc='CN SND: CARD CARD_1 ID C_0 C_1 C_2 C_3',
        precedence=11,
    )
    clay = param.DataFrame(
        default=pd.DataFrame({'ID': [], 'C_0': [], 'C_1': [], 'C_2': [], 'C_3': []}),
        doc='CN CLA: CARD CARD_1 ID C_0 C_1 C_2 C_3',
        precedence=12,
    )
    # TODO RCT, noqa an sand and clay

    def __init__(self):
        super(ConstituentProperties, self).__init__()
        self._update_salinity()
        self._update_temperature()
        self._update_vorticity()

    @param.depends('salinity', watch=True)
    def _update_salinity(self):
        self.param.salinity_id.precedence = -1
        self.param.reference_concentration.precedence = -1
        if self.salinity:
            self.param.salinity_id.precedence = 1.1
            self.param.reference_concentration.precedence = 2

    @param.depends('temperature', watch=True)
    def _update_temperature(self):
        self.param.temperature_id.precedence = -1
        self.param.reference_temperature.precedence = -1
        self.param.air_water_heat_transfer.precedence = -1
        self.param.short_wave_radiation_series.precedence = -1
        self.param.dew_point_temperature_series.precedence = -1
        if self.temperature:
            self.param.temperature_id.precedence = 3.1
            self.param.reference_temperature.precedence = 4
            self.param.air_water_heat_transfer.precedence = 5
            self.param.short_wave_radiation_series.precedence = 5.1
            self.param.dew_point_temperature_series.precedence = 5.2

    @param.depends('vorticity', watch=True)
    def _update_vorticity(self):
        self.param.vorticity_id.precedence = -1
        self.param.vorticity_normalization.precedence = -1
        self.param.vorticity_as_term.precedence = -1
        self.param.vorticity_ds_term.precedence = -1
        if self.vorticity:
            self.param.vorticity_id.precedence = 6.1
            self.param.vorticity_normalization.precedence = 7
            self.param.vorticity_as_term.precedence = 8
            self.param.vorticity_ds_term.precedence = 9

    def set_not_required(self, value=False):
        self.salinity = value
        self.temperature = value
        self.air_water_heat_transfer = value
        self.vorticity = value

    def panel(self):
        return pn.Pane(self.param, show_name=False)
