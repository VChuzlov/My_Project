import numpy as np
from flow import flow
import copy

class valve:
    def __init__(self, name: str, inlet_flow: flow, pressure_drop) -> None:
        self.name = name
        self.outlet_flow = copy.deepcopy(inlet_flow)
        self.pressure_drop = pressure_drop
        return
    
    def calculate(self):
        self.outlet_flow.pressure -= self.pressure_drop
        #Здесь должно быть распределение фаз по потокам
        self.outlet_flow.calculate()
        return self.outlet_flow

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
