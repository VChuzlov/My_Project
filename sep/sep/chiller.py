import numpy as np
from flow import flow
import copy

class chiller:
    def __init__(self, name: str, inlet_flow: flow, key_value, key_value_type = "temperature_out") -> None:
        self.name = name
        self.inlet_flow = copy.deepcopy(inlet_flow)
        self.key_value = key_value
        self.key_value_type = key_value_type
        return
    
    def calculate(self):
        self.outlet_flow = copy.deepcopy(self.inlet_flow)
        if self.key_value_type == "temperature_out":
            self.outlet_flow.temperature = self.key_value
        else: 
            self.outlet_flow.temperature = np.random.randint(-50, 50)
        self.delta_T = self.outlet_flow.temperature - self.inlet_flow.temperature
        self.duty = np.random.randint(0, 1000)
        #Здесь должно быть распределение фаз по потокам
        self.outlet_flow.calculate()
        return self.outlet_flow, self.duty

    def results(self):
        res_list = ['delta_T', 'duty']
        res = {}
        for property in res_list:
            res[f'{self.name}_{property}'] = self.__dict__[property]
        return res

def main() -> None:
    
    return 


if __name__ == '__main__':
    main()
