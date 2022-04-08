"""
    The Description of Peng-Robinson model.
"""

import constants as const
import numpy as np
from scipy.optimize import fsolve, root_scalar


def celsius_to_rankin(t):
    return t * 1.8 + 491.67


def celsius_to_kelvin(t):
    return t + 273.15


def kpa_to_psi(p):
    return p * 0.145


def get_tr(t, tc):
    return [t / tci for tci in tc]


def get_pr(p, pc):
    return [p / pci for pci in pc]


def solve_cubic_equation(a, b, c):
    """
    x ** 3 + a * x ** 2 + b * x + 1 = 0
    :param a:
    :param b:
    :param c:
    :return:
    """

    def sgn(x):
        if x > 0:
            return 1
        if x == 0:
            return 0
        if x < 0:
            return -1

    x1, x2, x3 = [None] * 3

    q = (a ** 2 - 3 * b) / 9
    r = (2 * a ** 3 - 9 * a * b + 27 * c) / 54
    s = q ** 3 - r ** 2

    if s > 0:
        fi = 1 / 3 * np.arccos(r / q ** (3 / 2))
        x1 = -2 * q ** 0.5 * np.cos(fi) - a / 3
        x2 = -2 * q ** 0.5 * np.cos(fi + 2 / 3 * 3.14) - a / 3
        x3 = -2 * q ** 0.5 * np.cos(fi - 2 / 3 * 3.14) - a / 3

    elif s < 0:
        if q > 0:
            fi = 1 / 3 * np.arccosh(abs(r) / q ** (3 / 2))
            x1 = -2 * sgn(r) * q ** 0.5 * np.cosh(fi) - a / 3
        elif q < 0:
            sss = abs(r) / abs(q) ** (3 / 2)
            fi = 1 / 3 * np.arcsinh(sss)
            x1 = -2 * sgn(r) * abs(q) ** 0.5 * np.sinh(fi) - a / 3
        else:
            x1 = -(c - a * a * a / 27) ** (1 / 3) - a / 3

    else:
        x1 = -2 * sgn(r) * q ** 0.5 - a / 3
        x2 = sgn(r) * q ** 0.5 - a / 3

    x1 = None if not x1 or x1 < 0 else x1
    x2 = None if not x2 or x2 < 0 else x2
    x3 = None if not x3 or x3 < 0 else x3

    x = [x for x in [x1, x2, x3] if x]

    return x


def p_sat_by_wilson(t, pc, tc, omega):
    """
    Wilson equation for Psat calculation/
    :param t: K
    :param pc: Pa
    :param tc: K
    :param omega:
    :return:
    """
    return [pci * np.exp(5.372697 * (1 + omi) * (1 - tci / t))
            for pci, omi, tci in zip(pc, omega, tc)]


def rachford_rice(e, zi, ki):
    """
    Rachford-Rice equation
    :param e: mole fraction of vapor phase
    :param zi: components mole fractions
    :param ki: phase equilibrium constant
    :return:
    """
    return sum(z * (k - 1) / (1 + e * (k - 1))
               for z, k in zip(zi, ki) if k > 0)


def first_guess(t, p, zi, pc, tc, omega):
    psat_i = p_sat_by_wilson(t, pc, tc, omega)

    ki = [psat / p for psat in psat_i]

    e0 = np.array([0.0])
    # e, *_ = fsolve(rachford_rice, e0, args=(zi, ki))
    e = root_scalar(rachford_rice, method='bisect', bracket=(0, 1), args=(zi, ki)).root
    e = 1 if e > 1 else e
    e = 0 if e < 0 else e

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

    return xi, yi, e, ki


def for_loop(zi, ki):
    e0 = np.array([0.0])
    e, *_ = fsolve(rachford_rice, e0, args=(zi, ki))
    e = 1 if e > 1 else e
    e = 0 if e < 0 else e

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

    return xi, yi, e


def second_guess(xi, yi, tr, pr, omega, kij):
    a0, b0 = 0.457235529, 0.077796

    bi = [b0 * pri / tri for pri, tri in zip(pr, tr)]
    bv = sum(y * b for y, b in zip(yi, bi))
    bl = sum(x * b for x, b in zip(xi, bi))

    ni = [0.37464 + 1.54226 * omi - 0.26992 * omi ** 2
          if omi <= 0.49
          else 0.379642 + (1.48503 - (0.164423 - 1.016666 * omi) * omi) * omi
          for omi in omega]

    alphai = [(1 + n * (1 - tri ** 0.5)) ** 2
              for n, tri in zip(ni, tr)]

    aci = [a0 * alph * pri / tri ** 2 for alph, pri, tri in zip(alphai, pr, tr)]

    # ai = [a * alph for a, alph in zip(aci, alphai)]
    qij = [
        [(1 - kij[i][j]) * (aci[i] * aci[j]) ** 0.5
         for j in range(const.COMP_COUNT)
         ]
        for i in range(const.COMP_COUNT)
    ]

    av, al = 0, 0
    for i in range(const.COMP_COUNT):
        for j in range(const.COMP_COUNT):
            av += yi[i] * yi[j] * qij[i][j]
            al += xi[i] * xi[j] * qij[i][j]

    return aci, bi, al, av, bl, bv, qij


def cubic_pr(z, a, b, c):
    return z ** 3 + a * z ** 2 + b * z + c


def get_z(foo, a, b, c, flag):
    z0 = np.array([0.1, 0.5, 1.0])
    sol = fsolve(foo, z0, args=(a, b, c))

    if flag == 'l':
        return sol[np.where(sol > 0)].min()

    if flag == 'v':
        return sol[np.where(sol > 0)].max()


def calculate_components_fugacity(mole_frac, bi, b, a, z, qij):
    phi = [0 for _ in range(const.COMP_COUNT)]
    fi = phi[:]

    for i in range(const.COMP_COUNT):
        s = 0

        for j in range(const.COMP_COUNT):
            s += mole_frac[j] * qij[i][j]

        if b:
            phi[i] = np.exp((z - 1) * bi[i] / b - np.log(z - b) -
                            a / (2 ** 0.5 * b) * (s / a - bi[i] / b / 2)
                            * np.log((z + (1 + 2 ** 0.5) * b) / (z - (-1 + 2 ** 0.5) * b)))

        fi[i] = phi[i] * mole_frac[i]

    return phi


def condition(zi, yi, xi, e, ki_prev, ki):
    cond2 = sum(z - y * e - (1 - e) * x for z, y, x in zip(zi, yi, xi))
    cond3 = round(sum(xi), 4) == 1
    cond4 = round(sum(yi), 4) == 1
    cond5 = sum(abs(k_prev - k) for k_prev, k in zip(ki_prev, ki)) <= 1e-4

    return cond5 and cond2 and cond3 and cond4


def calculate_equilibrium_by_pr(zi, t, p, tc, pc, omega, kij, foo=cubic_pr, cond=condition):
    t_r = celsius_to_rankin(t)
    tcr = [celsius_to_rankin(t_) for t_ in tc]

    p_psi = kpa_to_psi(p)
    pc_psi = [kpa_to_psi(p_) for p_ in pc]

    tr = get_tr(t_r, tcr)
    pr = get_pr(p_psi, pc_psi)

    xi, yi, e, ki_prev = first_guess(
        t_r, p_psi, zi, pc_psi, tcr, omega)
    i = 0
    while True:
        aci, bi, al, av, bl, bv, qij = second_guess(
            xi, yi, tr, pr, omega, kij
        )

        # zv = get_z(
        #     foo, flag='v',
        #     a=(bv - 1),
        #     b=(av - 2 * bv - 3 * bv ** 2),
        #     c=((-av + bv ** 2 + bv) * bv)
        # )
        zv, *_ = solve_cubic_equation(
            a=(bv - 1),
            b=(av - 2 * bv - 3 * bv ** 2),
            c=((-av + bv ** 2 + bv) * bv)
        )
        # zl = get_z(
        #     foo, flag='l',
        #     a=(bl - 1),
        #     b=(al - 2 * bl - 3 * bl ** 2),
        #     c=((-al + bl ** 2 + bl) * bl)
        # )
        zl, *_ = solve_cubic_equation(
            a=(bl - 1),
            b=(al - 2 * bl - 3 * bl ** 2),
            c=((-al + bl ** 2 + bl) * bl)
        )

        phiv = calculate_components_fugacity(
            yi, bi, bv, av, zv, qij
        )
        phil = calculate_components_fugacity(
            xi, bi, bl, al, zl, qij
        )

        ki = [pl / pv for pv, pl in zip(phiv, phil)]

        xi, yi, e, = for_loop(zi, ki)
        print(i, e)
        if cond(zi, yi, xi, e, ki_prev, ki):
            return xi, yi, e

        i += 1
        if i > 10000:
            return xi, yi, e

        ki_prev = ki[:]


class PRSolution:
    pass


if __name__ == '__main__':
    import flow

    # [2 * i for i in range(const.COMP_COUNT)]
    f = flow.Flow(mass_flows=const.mass_flows,
                  temperature=39.99,
                  pressure=12000)

    xi, yi, e = calculate_equilibrium_by_pr(
        f.mole_fractions, f.temperature,
        f.pressure,
        [tc for tc in const.TC],
        [pc for pc in const.PC],
        const.OMEGA, const.PR_Kij
    )

    print(e)
    print(xi)
