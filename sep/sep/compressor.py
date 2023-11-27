import numpy as np
from flow import flow
import copy

class compressor:
    def __init__(self, name: str, inlet_flow: flow,
                 key_value, key_value_type = "pressure_out",
                  efficiency_value = 75, efficiency_type = 'adiabatic') -> None:
        self.name = name
        self.inlet_flow = inlet_flow
        self.key_value = key_value
        self.key_value_type = key_value_type
        self.efficiency_value = efficiency_value
        self.efficiency_type = efficiency_type
        return
    
    def calculate(self):
        self.outlet_flow = copy.deepcopy(self.inlet_flow)
        if self.key_value_type == "pressure_out":
            self.outlet_flow.pressure = self.key_value
            self.delta_pressure = self.key_value - self.inlet_flow.pressure
        else: 
            self.delta_pressure = np.random.randint(0, 10000)
            self.outlet_flow.pressure = np.random.randint(self.inlet_flow.pressure, 10000)
        self.outlet_flow.calculate()
        self.compress_ratio = self.outlet_flow.pressure / self.inlet_flow.pressure
        self.duty = np.random.randint(0, 10000)

        self.head_poly = np.random.randint(0, 10000)
        self.head_adi = np.random.randint(0, 10000)
        self.efficiency_poly = np.random.randint(0, 100)
        self.efficiency_adi = self.efficiency_value
        self.polytropic_exp = np.random.randint(0, 100)
        self.adiabatic_exp =np.random.randint(0, 100)
        return self.outlet_flow, self.duty

    def results(self):
        res_list = ['key_value', 'key_value_type', 'efficiency_value', 'efficiency_type',
                    'delta_pressure', 'compress_ratio', 'duty', 'head_poly', 'head_adi',
                    'efficiency_poly', 'efficiency_adi', 'polytropic_exp', 'adiabatic_exp']
        res = {}
        for property in res_list:
            res[f'{self.name}_{property}'] = self.__dict__[property]
        return res

def main() -> None:
    
    return 


if __name__ == '__main__':
    main()
