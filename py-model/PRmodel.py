"""
    The Decribtion of Peng-Robinson model.
"""

import constants as const
import numpy as np
from scipy.optimize import fsolve


def get_component_pressure_by_pr(t, v, tc, pc, omega, kij, x):
    """
    Peng-Robinson equation for component pressure
    :param t: <float> temperature
    :param v: <float> molar volume
    :param tc: <float> critical temperature
    :param pc: <float> critical pressure
    :param omega: <float> acentric factor
    :return: <float> component pressure
    """
    b_i = [0.077796 * const.R * t / p for t, p in zip(tc, pc)]
    b = sum(b * xi for b, xi in zip(b_i, x))

    ac_i = [0.457235 * (const.R * t) ** 2 / p for t, p in zip(tc, pc)]
    m_i = [0.37646 + 1.54226 * om - 0.26992 * om ** 2 if om <= 0.49
           else 0.379642 + (1.48503 - (0.164423 - 1.016666 * om) * om) * om
           for om in omega]
    alpha_i = [(1 + m * (1 - (t / t_c) ** 0.5)) ** 2 for m, t_c in zip(m_i, tc)]
    a_i = [ac * alpha for ac, alpha in zip(ac_i, alpha_i)]

    aij = [[(1 - kij[i][j]) * (a_i[i] * a_i[j]) ** 0.5 for j in range(const.COMP_COUNT)]
           for i in range(const.COMP_COUNT)]

    a = 0
    for i in range(const.COMP_COUNT):
        for j in range(const.COMP_COUNT):
            a += aij[i][j] * x[i] * x[j]

    return [const.R * t / (vi - b) - a / (vi * (vi + b) + b * (vi - b))
            for vi in v]


def get_zfactor(t, v, tc, pc, omega, kij, x):
    """
    Peng-Robinson equation for Z-factor calculation
    :param t: <float> temperature
    :param v: <ndarray> molar volume
    :param tc: <ndarray> critical temperature
    :param pc: <ndarray> critical pressure
    :param omega: <ndarray> acentric factor
    :param kij: <ndarray> Peng-Robinson params
    :param x: <ndarray> components molar fractions
    :return: <float> Z-factor
    """
    b_i = 0.0077796 * const.R * tc / pc
    ac_i = 0.457235 * (const.R * tc) ** 2 / pc
    m_i = 0.37646 + 1.54226 * omega - 0.26992 * omega ** 2 if omega <= 0.49 \
        else 0.379642 + (1.48503 - (0.164423 - 1.016666 * omega) * omega) * omega
    alpha_i = 1 + m_i * (1 - (t / tс) ** 0.5) ** 2
    a_i = ac_i * alpha_i

    p = const.R * t / (v - b_i) - a_i / (v * (v + b_i) + b_i * (v - b_i))

    b = (x * b_i).sum()
    a = 0
    for i in range(x.shape[0]):
        s = 0
        for j in range(x.shape[0]):
            s += x[i] * x[j] + (a_i[i] * a_i[j]) ** 0.5 * (1 - kij[i][j])
        a += s

    A = a * p / (const.R * t) ** 2
    B = b * p / (const.R * t)

    equation = lambda z: z ** 3 + (1 - B) * z ** 2 \
                         + (A - 2 * B - 3 * B ** 2) * z \
                         - (A * B - B ** 2 - B ** 3)

    z = fsolve(equation, np.array([0, 0.5, 1]))
    return z.min()


class PRSolution:
    pass


if __name__ == '__main__':
    import flow

    f = flow.Flow(mass_flows=[2 * i for i in range(const.COMP_COUNT)],
                  temperature=273.15,
                  pressure=101.325 * 1e3)

    pi = get_component_pressure_by_pr(
            f.temperature,
            const.MOLAR_VOLUME,
            [tc + 273.15 for tc in const.TC],
            [pc * 1e3 for pc in const.PC],
            const.OMEGA,
            const.PR_Kij,
            f.mole_fractions
        )
    ki = [p / f.pressure for p in pi]

    s1 = sum(zf * k for zf, k in zip(f.mole_fractions, ki))
    s2 = sum(zf / k for zf, k in zip(f.mole_fractions, ki))

    xi = [0 for _ in range(const.COMP_COUNT)]
    yi = [0 for _ in range(const.COMP_COUNT)]
    e = 0

    foo = lambda x: sum(zf * (k - 1) / (1 + x * (k - 1))
                        for zf, k in zip(f.mole_fractions, ki))

    if s1 < 1:
        xi = f.mole_fractions[:]
    elif s2 < 1:
        yi = f.mole_fractions[:]
    else:
        e, *_ = fsolve(foo, np.array([0.5]))
        xi = [zf / (1 + e * (k - 1)) for zf, k in zip(f.mole_fractions, ki)]
        yi = [zf * k / (1 + e * (k - 1)) for zf, k in zip(f.mole_fractions, ki)]

    print(xi)
    print(sum(xi))
    print(yi)
    print(sum(yi))
    print(e)
