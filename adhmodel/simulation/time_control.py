# -*- coding: utf-8 -*-
import param
import panel as pn


class TimeControl(param.Parameterized):
    start_time = param.Number(
        default=0.0,
        bounds=(0, None),
        doc='TC T0: Starting time of the simulation. REQUIRED.',
        precedence=1,
    )
    start_time_units = param.ObjectSelector(
        default='None',
        objects=['None', '0 - seconds', '1 - minutes', '2 - hours', '3 - days', '4 - weeks'],
        doc='TC T0: Time units.',
        precedence=2,
    )
    end_time = param.Number(
        default=0.0,
        bounds=(0, None),
        doc='TC TF: Ending time of the simulation. REQUIRED.',
        precedence=3,
    )
    end_time_units = param.ObjectSelector(
        default='None',
        objects=['None', '0 - seconds', '1 - minutes', '2 - hours', '3 - days', '4 - weeks'],
        doc='TC TF: Time units.',
        precedence=4,
    )
    time_step_option = param.ObjectSelector(
        default='Steady state solution (TC STD)',
        objects=['Steady state solution (TC STD)', 'Time step series (SERIES DT)', 'Auto Time Step Find (TC ATF)'],
        doc="Time stepping option.",
        precedence=5,
    )
    max_time_step_size_time_series = param.Integer(
        default=1,
        bounds=(1, None),
        doc="SERIES DT: Series ID containing the series (time, max_time_step).",
        precedence=6,
    )
    steady_state_min_time_step_size = param.Number(
        default=0.0,
        bounds=(0, None),
        doc='TC STD: Minimum time step size.',
        precedence=7,
    )
    steady_state_max_time_step_size = param.Number(
        default=0.0,
        bounds=(0, None),
        doc='TC STD: Maximum time step size.',
        precedence=8,
    )
    auto_time_step_find_min_time_step_size = param.Number(
        default=0.0,
        bounds=(0, None),
        doc='TC ATF: Minimum time step size.',
        precedence=9,
    )
    auto_time_step_find_max_time_step_size_series = param.Integer(
        default=0,
        bounds=(0, None),
        doc='TC ATF: Maximum time step size series.',
        precedence=10,
    )
    #TODO TC NDP - need documentation

    def __init__(self):
        super(TimeControl, self).__init__()
        self._update_time_step_option()

    @param.depends('time_step_option', watch=True)
    def _update_time_step_option(self):
        self.param.max_time_step_size_time_series.precedence = -1
        self.param.steady_state_min_time_step_size.precedence = -1
        self.param.steady_state_max_time_step_size.precedence = -1
        self.param.auto_time_step_find_min_time_step_size.precedence = -1
        self.param.auto_time_step_find_max_time_step_size_series.precedence = -1
        if self.time_step_option == 'Steady state solution (TC STD)':
            self.param.steady_state_min_time_step_size.precedence = 7
            self.param.steady_state_max_time_step_size.precedence = 8
        elif self.time_step_option == 'Time step series (SERIES DT)':
            self.param.max_time_step_size_time_series.precedence = 6
        else:
            self.param.auto_time_step_find_min_time_step_size.precedence = 9
            self.param.auto_time_step_find_max_time_step_size_series.precedence = 10

    def panel(self):
        return pn.Pane(self.param, show_name=False)
