import numpy as np
from flow import flow
import copy

class separator2f:
    def __init__(self, name: str, inlet_flow: flow, pressure_drop = 0) -> None:
        self.name = name
        self.inlet_flow = copy.deepcopy(inlet_flow)
        self.pressure_drop = pressure_drop
        return
    
    def calculate(self):
        self.inlet_flow.pressure -= self.pressure_drop
        #Здесь должно быть распределение фаз по потокам
        self.outlet_gas_flow = copy.deepcopy(self.inlet_flow)
        self.outlet_liquid_flow = copy.deepcopy(self.inlet_flow)
        self.outlet_gas_flow.calculate()
        self.outlet_liquid_flow.calculate()
        return self.outlet_gas_flow, self.outlet_liquid_flow

    def results(self):
        res_list = ['pressure_drop']
        res = {}
        for property in res_list:
            res[f'{self.name}_{property}'] = self.__dict__[property]
        return res

def main() -> None:
    
    return 


if __name__ == '__main__':
    main()
