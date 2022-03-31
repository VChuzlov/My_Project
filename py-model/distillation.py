"""
    This module contains the description for DistillationColumn class.
"""
import numpy as np

import PRmodel as pr
from scipy.optimize import root_scalar


def get_t_sat(t, p, tc, pc, omega, v, model=pr.get_component_pressure_by_pr):
    ps = model(t, v, tc, pc, omega)
    ki = ps / p
    return ki - 1


class DistillationColumn:
    def __init__(self):
        self.initial_temperature_profile = []

    def __repr__(self):
        return

    def initial_guess_for_t_profile(self, fprofile, zfprofile, pressure_profile,
                                    tc, pc, omega, v, method='secant', method_args={}):
        p = np.median(pressure_profile)
        t_sat = [root_scalar(get_t_sat, method=method,
                             args=(p, t_c, p_c, omega_, v_, ),
                             **method_args)
                 for t_c, p_c, omega_, v_ in zip(tc, pc, omega, v)]

        sum_zf_times_f = [0 for _ in range(zfprofile.shape[0])]
        for i in range(zfprofile.shape[0]):
            sum_zf_times_f[i] = sum(zfprofile[i][j] * f for j, f in enumerate(fprofile))
        s = sum(t * s for t, s in zip(t_sat, sum_zf_times_f))

        sum_f = sum(f for f in fprofile)
        t_ave = s / sum_f
        for i in range(zfprofile.shape[0]):
            s +=
        self.initial_temperature_profile = t_sat

        return


if __name__ == '__main__':
    pass
