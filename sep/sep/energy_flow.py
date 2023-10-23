import numpy as np

class energy_flow:
    def __init__(self, name: str) -> None:
        self.name = name
        self.energy = None
        return
    

    def results(self):
        res_list = ['energy']
        res = {}
        for property in res_list:
            res[f'{self.name}_{property}'] = self.__dict__[property]
        return res
