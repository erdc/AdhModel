# -*- coding: utf-8 -*-
import param
import panel as pn


class IterationParameters(param.Parameterized):
    non_linear_iterations = param.Integer(
        default=10,
        bounds=(1, None),
        softbounds=(1, 20),
        doc='IP NIT: Maximum number of nonlinear iterations. REQUIRED.',
        precedence=1,
    )
    non_linear_tolerance_option = param.ObjectSelector(
        default='Specify residual and incremental (IP NTL & IP ITL)',
        objects=['Specify residual and incremental (IP NTL & IP ITL)', 'Specify residual (IP NTL)',
                 'Specify incremental (IP ITL)'],
        doc='Types of nonlinear tolerance(s) specified by user.',
        precedence=2,
    )
    non_linear_residual_tolerance = param.Number(
        default=0.001,
        bounds=(1e-20, None),
        softbounds=(0, 1),
        doc='IP NTL: Nonlinear residual tolerance. IP NTL and/or IP ITL REQUIRED.',
        precedence=3,
    )
    non_linear_incremental_tolerance = param.Number(
        default=0.01,
        bounds=(1e-20, None),
        softbounds=(0, 1),
        doc='IP ITL: Nonlinear incremental tolerance. IP NTL and/or IP ITL REQUIRED.',
        precedence=4,
    )
    linear_iterations = param.Integer(
        default=80,
        bounds=(1, None),
        softbounds=(1, 100),
        doc='IP MIT: Maximum number of linear iterations. REQUIRED.',
        precedence=5,
    )

    def __init__(self):
        super(IterationParameters, self).__init__()
        self._update_non_linear_residual_tolerance()

    @param.depends('non_linear_tolerance_option', watch=True)
    def _update_non_linear_residual_tolerance(self):
        self.param.non_linear_residual_tolerance.precedence = -1            # NTL
        self.param.non_linear_incremental_tolerance.precedence = -1         # ITL
        objects = list(self.param.non_linear_tolerance_option.get_range())
        if self.non_linear_tolerance_option == objects[0]:                      # both
            self.param.non_linear_residual_tolerance.precedence = 3
            self.param.non_linear_incremental_tolerance.precedence = 4
        elif self.non_linear_tolerance_option == objects[1]:                    # NTL only
            self.param.non_linear_residual_tolerance.precedence = 3
        else:                                                                   # ITL only
            self.param.non_linear_incremental_tolerance.precedence = 4

    def panel(self):
        return pn.Pane(self.param, show_name=False)
