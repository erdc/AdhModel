# -*- coding: utf-8 -*-
import param
import panel as pn


class MaterialWindProperties(param.Parameterized):
    stress_formulation = param.ObjectSelector(
        default=0,
        objects={'0 - no transform': 0, '1 - Wu': 1, '2 - Teeter': 2},
        doc='MP WND STR: Wind stress formulation.',
        precedence=1,
    )
    attenuation = param.Number(
        default=1.0,
        bounds=(0, 1),
        doc='MP WND ATT: Wind attenuation factor applied to wind speeds or shears.',
        precedence=2,
    )

    def panel(self):
        return pn.Pane(self.param, show_name=False)
