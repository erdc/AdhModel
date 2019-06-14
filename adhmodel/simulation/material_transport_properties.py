# -*- coding: utf-8 -*-
import param
import panel as pn


class MaterialTransportProperties(param.Parameterized):
    refinement_tolerance = param.Number(
        default=30.0,
        bounds=(0, None),
        doc='MP TRT: Transport refinement tolerance (constituent property).',
        precedence=2,
    )
    turbulent_diffusion_rate = param.Number(
        default=0.0,
        bounds=(0, None),
        doc='MP DF: Turbulent diffusion (constituent property).',
        precedence=3,
    )

    def panel(self):
        return pn.Pane(self.param, show_name=False)
