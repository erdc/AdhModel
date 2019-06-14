# -*- coding: utf-8 -*-
import param
import panel as pn


class OutputControl(param.Parameterized):
    output_control_option = param.ObjectSelector(
        default='Specify output frequency (OC)',
        objects=['Specify output frequency (OC)', 'Specify autobuild (SERIES AWRITE)'],
        doc='Output control option to specify output time using a time series (OC) or an autobuild time series (OS).',
        precedence=1,
    )
    oc_time_series_id = param.Integer(
        default=0,
        bounds=(0, None),
        doc='OC: Time Series ID for output times.',
        precedence=2,
    )
    output_flow_strings = param.DataFrame(
        default=None,
        doc='FLX: CARD S_ID, '
            'CARD - FLX, '
            'S_ID - String ID for mid string or edge string for which flow is to be output.',
    )
    print_adaptive_mesh = param.Boolean(
        default=False,
        doc='PC ADP: Adaptive mesh printing active.',
        precedence=4,
    )
    print_numerical_fish_surrogate = param.Boolean(
        default=False,
        doc='PC ELM: Print numerical fish surrogate information in TecPlot format.',
        precedence=5,
    )
    screen_output_residual = param.Boolean(
        default=False,
        doc='SOUT RESID: output the residual to the screen.',
        precedence=6,
    )
    screen_output_all = param.Boolean(
        default=False,
        doc='SOUT ALL: output all information to the screen.',
        precedence=7,
    )
    screen_output_mass_error = param.Boolean(
        default=False,
        doc='SOUT MERROR: Screen output mass error active.',
        precedence=8,
    )
    screen_output_worst_nonlinear_node = param.Boolean(
        default=False,
        doc='SOUT NLNODE: output the id of the worst nonlinear node.',
        precedence=9,
    )
    screen_output_worst_linear_node = param.Boolean(
        default=False,
        doc='SOUT LNODE: output the id of the worst linear node.',
        precedence=10,
    )
    file_output_wind = param.Boolean(
        default=False,
        doc='FOUT WIND: output wind to a file.',
        precedence=11,
    )
    file_output_wave = param.Boolean(
        default=False,
        doc='FOUT WAVE: output wave to a file.',
        precedence=12,
    )
    file_output_adapted_grid = param.Boolean(
        default=False,
        doc='FOUT ADAPT GRID: write adapted grid to a file.',
        precedence=13,
    )
    file_output_adapted_solution = param.Boolean(
        default=False,
        doc='FOUT ADAPT SW: write the adapted grid and solution to a file.',
        precedence=14,
    )
    file_output_adapted_transport = param.Boolean(
        default=False,
        doc='FOUT ADAPT CON: write the adapted transport files (does not include sediment).',
        precedence=15,
    )
    file_output_sediment = param.Boolean(
        default=False,
        doc='FOUT SED: write the adapted sediment files.',
        precedence=16,
    )

    def panel(self):
        return pn.Pane(self.param, show_name=False)
