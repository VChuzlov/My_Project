"""
    This module contains the description for DistillationColumn class.
"""
import numpy as np
import PRmodel as pr
from scipy.optimize import root_scalar
import json
import constants as const


def get_t_sat(t, p, tc, pc, omega, v, model=pr.get_component_pressure_by_pr):
    ps = model(t, v, tc, pc, omega)
    ki = ps / p
    return ki - 1


def calculation(rb, rd, wd, ld, ntrays):
    fj, uj, wj, lj0, vj0 = [], [], [], [], []

    # Lj0[1] := L0;
    lj0[1] = rd * ld
    for j in range(1, ntrays):
        lj0[j] = lj0[j - 1] + fj[j] - uj[j]
    lj0[ntrays - 1] = 9.3  # нужен пересчет

    for j in range(1, ntrays):
        vj0[j] = wd + ld + lj0[1] - wj[j]

    dj = [0 for _ in range(ntrays)]
    qj = [1 for _ in range(ntrays)]  # по идее надо его считать через энтальпию в зависимсоти от состояния
    for j in range(ntrays):
        dj[j] = (uj[j] + wj[j]) / (vj0[j] + lj0[j])

    vj0[ntrays - 1] = rb * lj0[ntrays - 1]
    for j in range(ntrays - 1, 1, -1):
        vj0[j] = ((1 - qj[j]) * fj[j] + vj0[j + 1]) / (dj[j] + 1)

    for j in range(1, ntrays - 1):
        lj0[j] = (fj[j] + vj0[j + 1] + lj0[j - 1]) / (dj[j] + 1) - vj0[j]

    s = 0
    for j in range(1, ntrays - 1):
        s = s + (qj[j] + rd) * fj[j] + (rd + 1) * uj[j] - wj[j]
    b = (rd * fj[1] + (rd + 1) * fj[ntrays - 1] + s) / (rd + rb + 1)  # {9.3}

    s = 0
    for j in range(1, ntrays - 1):
        s = s + (rb + 1 - qj[j]) * fj[j] + (rd + 1) * uj[j] + wj[j]
    d = ((rb + 1) * fj[1] + rd * fj[ntrays - 1] + s) / (rd + rb + 1)  # {4.5}

    return {'lj0': lj0, 'vj0': vj0, 'uj': uj, 'wj': wj, 'b': b, 'd': d}


def get_tray_material_balance_error(fj, uj, wj, lj, vj, ntrays):
    input_ = [0 for _ in range(ntrays)]
    output = [0 for _ in range(ntrays)]

    input_[1] = fj[1] + vj[2]
    input_[ntrays-1] = fj[ntrays-1] + lj[ntrays-2]
    output[ntrays-1] = wj[ntrays-1] + lj[ntrays-1] + vj[ntrays-1]

    for j in range(1, ntrays-1):
        input_[j] = fj[j] + lj[j-1] + vj[j+1]

    for j in range(ntrays-1):
        output[j] = uj[j] + wj[j] + lj[j] + vj[j]

    tray_errors = [inp - outp for inp, outp in zip(input_, output)]
    return tray_errors


def get_rb(a, b, fj, uj, wj, lj0,
           vj0, rd, wd, ld, ntrays, eps=1e-4):
    def f(rb):
        res = calculation(rb, rd, wd, ld, ntrays)
        args = [res[key] for key in ['uj', 'wj', 'lj0', 'vj0']]
        return get_tray_material_balance_error(fj, *args, ntrays)[0]

    if f(a) * f(b) < 0:
        while True:
            rb = (a + b) / 2
            if f(a) * f(rb) < 0:
                b = rb
            else:
                a = rb
            if (f(rb) == 0) or (abs(a - b) <= eps):
                break
    else:
        print('There is not solutions for rB!')

    return rb


class DistillationColumn:
    def __init__(self, trays_number, feed_tray):
        self.trays_number = trays_number
        self.feed_tray = feed_tray

        self.initial_temperature_profile = []
        self.lj0 = []
        self.vj0 = []

    def __repr__(self):
        column = {
            'Number of trays': self.trays_number,
            'Feed Tray': self.feed_tray,
            'Tj0': self.initial_temperature_profile,
            'Lj0': self.lj0,
            'Vj0': self.vj0,
        }
        json_string = json.dumps(column, indent=4)
        return json_string

    def t_profile_initial_guess_for(self, fprofile, zfprofile, pressure_profile,
                                    tc, pc, omega, v, method='secant', method_args={}):
        p = np.median(pressure_profile)
        t_sat = [root_scalar(get_t_sat, method=method, x0=10, x1=1000,
                             args=(p, t_c, p_c, omega_, v_,),
                             **method_args).root
                 for t_c, p_c, omega_, v_ in zip(tc, pc, omega, v)]

        sum_zf_times_f = [0 for _ in range(zfprofile.shape[0])]
        for i in range(zfprofile.shape[0]):
            sum_zf_times_f[i] = sum(zfprofile[i][j] * f for j, f in enumerate(fprofile))
        s = sum(t * s for t, s in zip(t_sat, sum_zf_times_f))

        sum_f = sum(f for f in fprofile)
        t_ave = s / sum_f

        t_min = sum(abs(ts - t_ave) * s
                    for ts, s in zip(t_sat, sum_zf_times_f)) / sum_f

        t_0 = [t_min + 2 * j / zfprofile.shape[1] * (t_ave - t_min)
               for j in range(zfprofile.shape[1])]

        self.initial_temperature_profile = t_0

        return

    def lv_initial_guess(self, fprofile, wd, ld, l0):
        lj0 = [0 for _ in range(self.trays_number)]
        vj0 = [0 for _ in range(self.trays_number)]
        uj = [0 for _ in range(self.trays_number)]
        wj = [0 for _ in range(self.trays_number)]

        rd = l0 / ld
        rb = 3
        uj[0] = ld
        wj[0] = wd
        uj[self.trays_number-1] = fprofile[self.feed_tray-1] - (ld + wd)  # д.б.сумма Fj
        vj0[0] = 0  # должно зависеть от типа конденсатора
        lj0[0] = ld * rd

        res = calculation(rb, rd, wd, ld, self.trays_number)
        args = [res[key] for key in ['uj', 'wj', 'lj0', 'vj0']]
        rb = get_rb(1e-5, 1000, fprofile, *args, rd, wd, ld, self.trays_number)
        res = calculation(rb, rd, wd, ld, self.trays_number)

        self.lj0 = res['lj0']
        self.vj0 = res['vj0']

        return


if __name__ == '__main__':
    column = DistillationColumn(12, 5)

    ntrays = 12
    feed_tray = 5
    fprofile = [100 if i == feed_tray else 0 for i in range(ntrays)]
    z = [[1 for i in range(ntrays)] for _ in range(const.COMP_COUNT)]
    pressure_profile = [0.1 for _ in range(ntrays)]
    column.t_profile_initial_guess_for(
        fprofile, np.array(z), pressure_profile,
        const.TC, const.PC, const.OMEGA, const.V
    )
    print(column.initial_temperature_profile)
