"""
    The Description of Peng-Robinson model.
"""

import constants as const
import numpy as np
from scipy.optimize import fsolve


def p_sat_by_wilson(t, pc, tc, omega):
    """
    Wilson equation for Psat calculation/
    :param t: K
    :param pc: Pa
    :param tc: K
    :param omega:
    :return:
    """
    return [pci * np.exp(5.373 * (1 - omi) * (1 - tci / t))
            for pci, omi, tci in zip(pc, omega, tc)]


def rachford_rice(e, zi, ki):
    """
    Rachford-Rice equation
    :param e: mole fraction of vapor phase
    :param zi: components mole fractions
    :param ki: phase equilibrium constant
    :return:
    """
    return sum(z * (k - 1) / (e * (k - 1) + 1)
               for z, k in zip(zi, ki))


def first_guess(t, p, zi, pc, tc, omega):
    psat_i = p_sat_by_wilson(t, pc, tc, omega)

    ki = [psat / p for psat in psat_i]

    e0 = np.array([0.5])
    e, *_ = fsolve(rachford_rice, e0, args=(zi, ki))

    s1 = sum(z * k for z, k in zip(zi, ki))
    s2 = sum(z / k for z, k in zip(zi, ki))

    xi = [0 for _ in range(const.COMP_COUNT)]
    yi = [0 for _ in range(const.COMP_COUNT)]

    if s1 < 1:
        xi = zi[:]

    elif s2 < 1:
        yi = zi[:]

    else:
        xi = [z / (1 + e * (k - 1)) for z, k in zip(zi, ki)]
        yi = [z * k / (1 + e * (k - 1)) for z, k in zip(zi, ki)]

    return xi, yi


def second_guess(t, p, xi, yi, tc, pc, omega, kij):
    a0, b0 = 0.457235, 0.077796

    bi = [b0 * const.R * tci / pci for tci, pci in zip(tc, pc)]
    bv = sum(y * b for y, b in zip(yi, bi))
    bl = sum(x * b for x, b in zip(xi, bi))

    ni = [0.37646 + 1.5422 * omi - 0.26992 * omi ** 2 if omi <= 0.49
          else 0.379642 + (1.48503 - (0.164423 - 1.016666 * omi) * omi) * omi
          for omi in omega]

    alphai = [(1 + n * (1 - (t / tci) ** 0.5)) ** 2
              for n, tci in zip(ni, tc)]

    aci = [a0 * (const.R * tci) ** 2 / pci
           for tci, pci in zip(tc, pc)]

    ai = [a * alph for a, alph in zip(aci, alphai)]

    av, al = 0, 0
    for i in range(const.COMP_COUNT):
        for j in range(const.COMP_COUNT):
            av += yi[i] * yi[j] * (ai[i] * ai[j]) ** 0.5 * (1 - kij[i][j])
            al += xi[i] * xi[j] * (ai[i] * ai[j]) ** 0.5 * (1 - kij[i][j])

    aav = av * p * (const.R * t) ** 2
    aal = al * p / (const.R * t) ** 2

    bbv = bv * p / (const.R * t)
    bbvl = bl * p / (const.R * t)

    return aav, aal, bbv, bbvl


def cubic_pr(z, a, b):
    return (z ** 3 + (b - 1) * z ** 2
            + (a - 2 * b - 3 * b ** 2) * z
            + (b ** 2 + b ** 3 - a * b))


def get_z(foo, z, a, b, flag):
    z0 = np.array([0.1, 0.5, 1.0])
    sol = fsolve(foo, z0, args=(a, b))

    if flag == 'l':
        return sol[np.where(sol > 0)].min()

    if flag == 'v':
        return sol[np.where(sol > 0)].max()


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
