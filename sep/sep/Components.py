
import sys
sys.path.insert(0, '')

import SQLite_sub as sql
import json

class Component:

    '''
    Класс предназначен для хранения информации по компоненту
    '''
    def __init__(self, prop_comp, use_db = True, use_json = False):
        '''
        Инициализация компонента
        :param prop_comp: структура-словарь, в которой содержится основная информация о компоненте
        :param use_db: Флаг,
        '''
        self.Name_comp = prop_comp["Name_comp"]
        self.p_cr = prop_comp["p_cr"] # кПа
        self.T_cr = prop_comp["T_cr"] # 0С
        self.w = prop_comp["w"]
        self.Mr = prop_comp["Mr"] # г/моль
        self.Tkip = prop_comp["Tkip"] # С
        if use_db:
            res = sql.get_prop_of_name(self.Name_comp)
            self.p_cr = float(res["CriticalPressureValue"].replace(',','.'))
            self.T_cr = float(res["CriticalTemperatureValue"].replace(',','.'))
            self.w = float(res["AcentricityValue"].replace(',','.'))
            self.Mr = float(res["MolecularWeightValue"].replace(',','.'))
            self.Tkip = float(res["NormalBoilingPointValue"].replace(',','.'))
            self.V_cr = float(res["CriticalVolumeValue"].replace(',','.'))
            self.ro_liq = float(res["StdLiquidDensityValue"].replace(',','.'))
            self.IH_A = float(res["IdealHCoeffsAValue"].replace(',', '.'))
            self.IH_B = float(res["IdealHCoeffsBValue"].replace(',', '.'))
            self.IH_C = float(res["IdealHCoeffsCValue"].replace(',', '.'))
            self.IH_D = float(res["IdealHCoeffsDValue"].replace(',', '.'))
            self.IH_E = float(res["IdealHCoeffsEValue"].replace(',', '.'))
            self.IH_F = float(res["IdealHCoeffsFValue"].replace(',', '.'))
            self.IH_G = float(res["IdealHCoeffsGValue"].replace(',', '.'))
            self.IH_H = float(res["IdealHCoeffsHValue"].replace(',', '.'))
            self.A_A = float(res["AntoineCoeffsAValue"].replace(',', '.'))
            self.A_B = float(res["AntoineCoeffsBValue"].replace(',', '.'))
            self.A_C = float(res["AntoineCoeffsCValue"].replace(',', '.'))
            self.A_D = float(res["AntoineCoeffsDValue"].replace(',', '.'))
            self.A_E = float(res["AntoineCoeffsEValue"].replace(',', '.'))
            self.A_F = float(res["AntoineCoeffsFValue"].replace(',', '.'))
            self.A_G = float(res["AntoineCoeffsGValue"].replace(',', '.'))
            self.CAS = str(res['CAS_Number2'].replace(',', '.'))
            self.names = str(res['TaggedName'].replace(',', '.'))
            self.H_obr = float(res['HeatOfFormationValue'].replace(',', '.'))


        if use_json:
            s_json = prop_comp
            self.Name_comp = s_json["Name_comp"]
            self.p_cr = float(s_json["p_cr"])
            self.T_cr = float(s_json["T_cr"])
            self.w = float(s_json["w"])
            self.Mr = float(s_json["Mr"])
            self.Tkip = float(s_json["Tkip"])
            self.V_cr = float(s_json["V_cr"])
            self.ro_liq = float(s_json["ro_liq"])
            self.IH_A = float(s_json["IH_A"])
            self.IH_B = float(s_json["IH_B"])
            self.IH_C = float(s_json["IH_C"])
            self.IH_D = float(s_json["IH_D"])
            self.IH_E = float(s_json["IH_E"])
            self.IH_F = float(s_json["IH_F"])
            self.IH_G = float(s_json["IH_G"])
            self.IH_H = float(s_json["IH_H"])
            self.p_cr_am = self.p_cr * 0.145
            self.T_cr_am = self.T_cr * 1.8 + 491.67


    def get_json(self):
        x = {
            "Name_comp": self.Name_comp,
            "p_cr": self.p_cr,
            "T_cr": self.T_cr,
            "w": self.w,
            "Mr": self.Mr,
            "Tkip": self.Tkip,
            "V_cr": self.V_cr,
            "ro_liq": self.ro_liq,
            "IH_A": self.IH_A,
            "IH_B": self.IH_B,
            "IH_C": self.IH_C,
            "IH_D": self.IH_D,
            "IH_E": self.IH_E,
            "IH_F": self.IH_F,
            "IH_G": self.IH_G,
            "IH_H": self.IH_H,
            "p_cr_am": self.p_cr_am,
            'T_cr_am': self.T_cr_am
        }
        # convert into JSON:
        y = json.dumps(x)

        # the result is a JSON string:
        return y

if __name__ == "__main__":
    com = Component({'Name_comp': 'Methane', 'p_cr': 4640.680176, 'T_cr': -82.45099487, 'w': 0.0114984, 'Mr': 16.04290009, 'Tkip': -161.525, 'V_cr': 0.098999903, 'ro_liq': 299.3940125, 'IH_A': -12.98, 'IH_B': 2.36459, 'IH_C': -0.00213247, 'IH_D': 5.66e-06, 'IH_E': -3.72e-09, 'IH_F': 8.61e-13, 'IH_G': 1.0, 'IH_H': 0.0, 'p_cr_am': 672.8986255199999, 'T_cr_am': 343.258209234},
                    False, True)
    print(com.IH_E)



