"""
    This module contains the description of class Flow
"""

import regularfunctions as rf
import constants as const
import json


class Flow:
    def __init__(self, mass_flow=None, std_liq_volume_flow=None, mole_flow=None,
                 mass_fractions=(), volume_fractions=(), mole_fractions=(),
                 mass_flows=(), volume_flows=(), mole_flows=(),
                 temperature=None, pressure=None,
                 liquid_density=None):
        if not (
                (mass_flow or mole_flow or std_liq_volume_flow) and
                (
                        any(mass_fractions) or any(mole_fractions) or any(volume_fractions) or
                        any(mass_flows) or any(mole_flows) or any(volume_flows)
                )
        ) and not temperature and not pressure:
            print('Set flow parameters')
            return

        # temperature and pressure
        self.temperature = temperature
        self.pressure = pressure

        # flow rates
        self.mass_flow = mass_flow
        self.std_liquid_volume_flow = std_liq_volume_flow
        self.mole_flow = mole_flow

        # flow composition
        if len(mass_fractions) == const.COMP_COUNT:
            self.mass_fractions = rf.norm(mass_fractions)
            self.mole_fractions = rf.convert_mass_to_mole_fractions(
                self.mass_fractions, const.MR)
            self.volume_fractions = rf.convert_mass_to_volume_fractions(
                self.mass_fractions, const.DENSITY)

        elif len(mole_fractions) == const.COMP_COUNT:
            self.mole_fractions = rf.norm(mole_fractions)
            self.mass_fractions = rf.convert_mole_to_mass_fractions(
                self.mole_fractions, const.MR)
            self.volume_fractions = rf.convert_mole_to_volume_fractions(
                self.mole_fractions, const.MR, const.DENSITY)

        elif len(volume_fractions) == const.COMP_COUNT:
            self.volume_fractions = rf.norm(volume_fractions)
            self.mass_fractions = rf.convert_volume_to_mass_fractions(
                self.volume_fractions, const.DENSITY)
            self.mole_fractions = rf.convert_volume_to_mole_fractions(
                self.volume_fractions, const.MR, const.DENSITY)

        elif len(mass_flows) == const.COMP_COUNT:
            self.mass_flows = mass_flows
            self.mass_fractions = rf.convert_flows_to_fractions(mass_flows)
            self.mass_flow = sum(self.mass_flows)
            self.mole_fractions = rf.convert_mass_to_mole_fractions(
                self.mass_fractions, const.MR)
            self.volume_fractions = rf.convert_mass_to_volume_fractions(
                self.mass_fractions, const.DENSITY)

        elif len(volume_flows) == const.COMP_COUNT:
            self.volume_flows = volume_flows
            self.volume_fractions = rf.convert_flows_to_fractions(volume_flows)
            self.std_liquid_volume_flow = sum(self.volume_flows)
            self.mass_fractions = rf.convert_volume_to_mass_fractions(
                self.volume_fractions, const.DENSITY)
            self.mole_fractions = rf.convert_volume_to_mole_fractions(
                self.volume_fractions, const.MR, const.DENSITY)

        elif len(mole_flows) == const.COMP_COUNT:
            self.mole_flows = mole_flows
            self.mole_fractions = rf.convert_flows_to_fractions(self.mole_flows)
            self.mole_flow = sum(self.mole_flows)
            self.mass_fractions = rf.convert_mole_to_mass_fractions(
                self.mole_fractions, const.MR)
            self.volume_fractions = rf.convert_mole_to_volume_fractions(
                self.mole_fractions, const.MR, const.DENSITY)

        # flow rates
        if self.mass_flow:
            self.volume_flows = rf.convert_mass_to_volume_flows(
                self.mass_flows, const.DENSITY)
            self.mole_flows = rf.convert_mass_to_mole_flows(
                self.mass_flows, const.MR)

            self.std_liquid_volume_flow = rf.get_volume_flow_rate(
                self.volume_flows)
            self.mole_flow = rf.get_mole_flow_rate(self.mole_flows)

        elif self.std_liquid_volume_flow:
            self.mass_flows = rf.convert_volume_to_mass_flows(
                self.volume_flows, const.DENSITY)
            self.mole_flows = rf.convert_volume_to_mole_flows(
                self.volume_flows, const.DENSITY, const.MR)

            self.mass_flow = rf.get_mass_flow_rate(mass_flows)
            self.mole_flows = rf.get_mole_flow_rate(self.mole_flows)

        elif self.mole_flow:
            self.mass_flows = rf.convert_mole_to_mass_flows(
                self.mole_flows, const.MR)
            self.volume_flows = rf.convert_mole_to_volume_flows(
                self.mole_flows, const.DENSITY, const.MR)

            self.mass_flow = rf.get_mass_flow_rate(self.mass_flows)
            self.std_liquid_volume_flow = rf.get_volume_flow_rate(self.volume_flows)

        # flow density
        if not liquid_density:
            self.liquid_density = rf.get_flow_density(
                self.mass_fractions, const.DENSITY)
        else:
            self.liquid_density = liquid_density

    def __repr__(self):
        flow = {
            'Temperature': self.temperature,
            'Pressure': self.pressure,

            'Mass Flow Rate': self.mass_flow,
            'Std Liquid Flow Rate': self.std_liquid_volume_flow,
            'Mole Flow Rate': self.mole_flow,

            'Mass Fractions': dict(zip(const.NAMES, self.mass_fractions)),
            'Volume Fractions': dict(zip(const.NAMES, self.volume_fractions)),
            'Mole Fractions': dict(zip(const.NAMES, self.mole_fractions)),
        }
        jsonstring = json.dumps(flow, indent=4)
        return jsonstring

    def mixing(self, flows_list):
        """
        Mixing flows from flows_list. Return resulting Flow object.
        :param flows_list: <list-object> A list of Flow object to mix.
        :return: Resulting Flow object.
        """
        mass_flow_list = [flow.mass_flow for flow in flows_list] + [self.mass_flow]
        mass_fractions_list = [flow.mass_fractions for flow in flows_list] + [self.mass_fractions]
        temperature_list = [flow.temperature for flow in flows_list] + self.temperature
        pressure_list = [flow.pressure for flow in flows_list] + self.pressure

        mixture = rf.mix_flows(mass_flow_list, mass_fractions_list,
                               temperature_list, pressure_list)
        return Flow(**mixture)


if __name__ == '__main__':
    f = Flow(
        mass_flows=[f * 10 for f in range(const.COMP_COUNT)],
        temperature=273.15,
        pressure=0.1,
        liquid_density=0.5
    )

    print(f)
