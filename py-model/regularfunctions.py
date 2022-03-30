import constants as const
from scipy.optimize import fsolve


# conversion of flow composition
def convert_mass_to_mole_fractions(mass_fractions, mr):
    s = sum(mf / m for mf, m in zip(mass_fractions, mr))
    return [mf / m / s for mf, m in zip(mass_fractions, mr)]


def convert_mass_to_volume_fractions(mass_fractions, density):
    s = sum(mf * d for mf, d in zip(mass_fractions, density))
    return [mf * d / s for mf, d in zip(mass_fractions, density)]


def convert_mole_to_mass_fractions(mole_fractions, mr):
    s = sum(m * mf for m, mf in zip(mr, mole_fractions))
    return [mf * m / s for mf, m in zip(mole_fractions, mr)]


def convert_mole_to_volume_fractions(mole_fractions, mr, density):
    s = sum(mf * m / d for mf, m, d in zip(mole_fractions, mr, density))
    return [mf * m / d / s for mf, m, d in zip(mole_fractions, mr, density)]


def convert_volume_to_mass_fractions(volume_fractions, density):
    s = sum(vf * d for vf, d in zip(volume_fractions, density))
    return [vf * d / s for vf, d in zip(volume_fractions, density)]


def convert_volume_to_mole_fractions(volume_fractions, mr, density):
    s = sum(vf * d / m for vf, d, m in zip(volume_fractions, density, mr))
    return [vf * d / m / s for vf, d, m in zip(volume_fractions, density, mr)]


# composition normalization
def norm(composition):
    s = sum(composition)

    return [x / s for x in composition]


# component flows
def get_mass_flows(mass_flow_rate, wt_fractions):
    return [mass_flow_rate * fraction for fraction in wt_fractions]


def get_mole_flows(mole_flow_rate, mole_fractions):
    return [mole_flow_rate * fraction for fraction in mole_fractions]


def get_volume_flows(volume_flow_rate, volume_fractions):
    return [volume_flow_rate * fraction for fraction in volume_fractions]


def convert_mole_to_mass_flows(mole_flows, mr):
    """
        mole_flows [kmole / h] * mr [g / mole = kg / kmole] = mass_flows [kg / h]
    """
    return [mf * m for mf, m in zip(mole_flows, mr)]


def convert_mole_to_volume_flows(mole_flows, density, mr):
    """
        mole_flows [kmole / h] / density [kg / m**3] * mr [g / mole = kg / kmole] = volume_flows [m**3 / h]
    """
    return [mf / (d * 1000) * m for mf, d, m in zip(mole_flows, density, mr)]


def convert_mass_to_volume_flows(mass_flows, density):
    """
        mass_flows [kg / h] / density [kg / m**3] = volume_flows [m**3 / h]
    """
    return [mf / (d * 1000) for mf, d in zip(mass_flows, density)]


def convert_mass_to_mole_flows(mass_flows, mr):
    """
        mass_flows [kg / h] / mr [g / mole = kg / kmole] = mole_flows [kmole / h]
    """
    return [mf / m for mf, m in zip(mass_flows, mr)]


def convert_volume_to_mass_flows(volume_flows, density):
    """
        volume_flows [m**3 / h] * density [kg / m**3] = mass_flows [kg / h]
    """
    return [vf * d * 1000 for vf, d in zip(volume_flows, density)]


def convert_volume_to_mole_flows(volume_flows, density, mr):
    """
        volume_flows [m**3 / h] * density [kg / m**3] / mr [g / mole = kg / kmole] = mole_flows [kmole / h]
    """
    return [vf * d * 1000 / m for vf, d, m in zip(volume_flows, density, mr)]


# flows to fractions conversion
def convert_flows_to_fractions(flows):
    flow_rate = sum(flows)
    return [flow / flow_rate for flow in flows]


# flow rates
def get_mole_flow_rate(mole_flows):
    return sum(mole_flows)


def get_mass_flow_rate(mass_flows):
    return sum(mass_flows)


def get_volume_flow_rate(volume_flows):
    return sum(volume_flows)


# flow average molecular weight
def get_flow_molecular_weight(mole_fractions, mr):
    return sum(m * mf for m, mf in zip(mr, mole_fractions))


# flow density
def get_flow_density(mass_fractions, density):
    s = sum(mf / d for mf, d in zip(mass_fractions, density))

    if not s == 0:
        return 1 / s

    return


# flow mass heat capacity
def get_mass_flow_cp(coef_cp, mass_fractions, temperature):
    """
        This function calculates flow heat capacity from mass heat capacities of individual components in dependency of
        temperature and flow mass composition. Coefficients is taken from enthalpy temperature dependency.
        Polynom produced by derivation.
    """
    k1 = [col[0] for col in coef_cp]
    k2 = [col[1] for col in coef_cp]
    k3 = [col[2] for col in coef_cp]
    k4 = [col[3] for col in coef_cp]
    k5 = [col[4] for col in coef_cp]

    def polynom(a, b, c, d, e, fraction):
        return (a + 2 * b * temperature + 3 * c * temperature ** 2 + 4 * d * temperature ** 3
                + 5 * e * temperature ** 4) * fraction

    result = sum(map(polynom, k1, k2, k3, k4, k5, mass_fractions))

    return result


def mix_flows(mass_flow_rates, mass_fractions, flow_temperatures, flow_presssures):
    """
        This function takes 3 lists with mass flow rates, mass compositions and flow temperatures.
        Returns a dict with mass_fractions, mass_flow_rate and flow temperature of the mixture.
    """
    results = {}

    # material balance

    flows_count = len(mass_flow_rates)
    components_count = len(mass_fractions[0])

    mixture_mass_flow_rate = sum(mass_flow_rates)

    components_mass_flows = [[mass_fractions[i][j] * mass_flow_rates[i] for i in range(flows_count)]
                             for j in range(components_count)]

    mixture_mass_fractions = [sum(components_mass_flows[i]) / mixture_mass_flow_rate for i in range(components_count)]

    # heat balance

    mass_flows_cp = [get_mass_flow_cp(const.COEF_MASS_CP, mass_fractions[i], flow_temperatures[i]) for i in
                     range(flows_count)]

    def equation(temperature):
        nonlocal mass_flows_cp, mass_flow_rates, flow_temperatures

        equation = sum(map(lambda a, b, c: a * b * c, mass_flows_cp, mass_flow_rates, flow_temperatures)) \
                   / (get_mass_flow_cp(const.COEF_MASS_CP, mixture_mass_fractions, temperature)
                      * mixture_mass_flow_rate) - temperature

        return equation

    initial_guess = sum(flow_temperatures) / flows_count
    mixture_temperature = fsolve(equation, initial_guess)

    results['mass_fractions'] = mixture_mass_fractions
    results['mass_flow'] = mixture_mass_flow_rate
    results['temperature'] = mixture_temperature[-1]
    results['pressure'] = min(flow_presssures)

    return results


if __name__ == '__main__':
    pass
