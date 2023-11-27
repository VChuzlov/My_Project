import numpy as np
from PengRobinson import PengRobinson

class flow:
    def __init__(self, name: str) -> None:
        self.name = name
        self.mass_flow = None
        self.temperature = None
        self.pressure = None
        self.composition = None
        self.package = None
        return
    
    def calculate(self):

        eos = self.package(composition=self.composition, pressure=self.pressure, temperature=self.temperature, mass_flow=self.mass_flow)
        eos.calculation()

        self.vapor_phase = eos.gas['Gm, kmole/h'] / (eos.gas['Gm, kmole/h'] + eos.liq['Gm, kmole/h'] + eos.water['Gm, kmole/h'])
        self.molar_flow = eos.gas['Gm, kmole/h'] + eos.liq['Gm, kmole/h'] + eos.water['Gm, kmole/h']
        self.volume_flow = self.molar_flow * 8.314 * (20+273) / self.pressure  # p*Gv = Gm*RT @20C + 1atm
        self.molar_enthalpy = eos.H
        self.mw_avg = eos.mw_avg
        self.density_mass = eos.ro
        self.z_factor = eos.Zv
        self.volume_flow_actual = self.mass_flow / self.density_mass
        self.Cp = eos.Cp
        self.Cv = eos.Cv
        self.Reid_pressure = eos.reid_pressure
        self.dew_point_water = eos.water_dew_point
        self.dew_point_hc = eos.hc_dew_point

        del eos
        return

    def results(self):
        res_list = ['mass_flow', 'temperature', 'pressure', 'vapor_phase',
                    'molar_flow', 'volume_flow', 'molar_enthalpy', 'mw_avg',
                    'density_mass', 'z_factor', 'volume_flow_actual',
                    'Cp', 'Cv', 'Reid_pressure', 'dew_point_water', 'dew_point_hc']
        res = {}
        for property in res_list:
            res[f'{self.name}_{property}'] = self.__dict__[property]
        return res

def main() -> None:
    flow1 = flow('flow1')
    flow1.pressure = 101.325
    flow1.temperature = -200
    flow1.mass_flow = 100
    flow1.composition = {
        'Methane': 1,
    }

    """flow1.composition = {
        'Methane': 0.5,
        'Ethane': 0.3,
        'Propane': 0.15,
        'H2O': 0.05
    }"""
    flow1.package = PengRobinson
    flow1.calculate()
    res = flow1.results()

    for property in res:
        print(f'{property}: {res[property]}')

    return 


if __name__ == '__main__':
    main()
