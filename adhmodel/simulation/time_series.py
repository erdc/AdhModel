# -*- coding: utf-8 -*-
import param
import panel as pn
import pandas as pd


class TimeSeries(param.Parameterized):
    series_id = param.Integer(
        default=1,
        bounds=(1, None),
        doc='DATASERIES: ID of this xy series.',
        precedence=1,
    )
    series_type = param.ObjectSelector(
        default='SERIES BC',
        objects=['SERIES AWRITE', 'SERIES WRITE', 'SERIES BC', 'SERIES DT', 'SERIES WIND', 'SERIES WAVE'],
        doc='Type of time series. Supported types in ADH 5.0 are "SERIES AWRITE", "SERIES WRITE", "SERIES BC", '
            '"SERIES DT", "SERIES WIND", "SERIES WAVE".',
        precedence=2,
    )
    time_series = param.DataFrame(
        default=pd.DataFrame(data=[], columns=['X', 'Y']),
        doc='Time series values.',
        precedence=3,
    )
    units = param.ObjectSelector(
        default=0,
        objects={'0 - seconds': 0, '1 - minutes': 1, '2 - hours': 2, '3 - days': 3, '4 - weeks': 4},
        doc='Time units.',
        precedence=4,
    )
    output_units = param.ObjectSelector(
        default=0,
        objects={'0 - seconds': 0, '1 - minutes': 1, '2 - hours': 2, '3 - days': 3, '4 - weeks': 4},
        doc='Output time units.',
        precedence=5,
    )
    x_location = param.Number(
        default=0,
        bounds=(None, None),
        doc="X location for SERIES WIND and SERIES WAVE.",
        precedence=6,
    )
    y_location = param.Number(
        default=0,
        bounds=(None, None),
        doc="Y location for SERIES WIND and SERIES WAVE.",
        precedence=7,
    )

    def __init__(self, **params):
        super(TimeSeries, self).__init__(**params)

    def panel(self):
        return pn.panel(self.param, show_name=False)
