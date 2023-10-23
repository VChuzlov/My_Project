import sys
sys.path.insert(0, '')
# import matplotlib.pyplot as plt
import math


class ComponentProperties:

    def __init__(self, comps) -> None:
        self.comps = comps
        all_comps = ['Methane', 'Ethane', 'Propane', 'i-Butane', 'n-Butane', 'i-Pentane', 'n-Pentane', 'n-Hexane', 'n-Heptane', 'n-Octane', 'H2O', 'Methanol']

        p_cr = [4640.680176, 4883.850098, 4256.660156, 3647.620117, 3796.620117, 3333.590088, 3375.120117, 3031.620117, 2736.780029, 2496.620117, 22120.0, 7376.450195]  # критические давления, кПа
        T_cr = [-82.45099487, 32.27800903, 96.74801025, 134.9460083, 152.0490051, 187.2480103, 196.4500061, 234.7480103, 267.00802, 295.4480225, 374.1490112, 239.4480225]  # критические температуры, С
        w = [0.0114984, 0.0986, 0.152400002, 0.18479, 0.201000005, 0.222240001, 0.253890008, 0.300700009, 0.349790007, 0.401800007, 0.344000012, 0.556990027]  # ацентрический фактор
        Mr = [16.04290009, 30.06990051, 44.09700012, 58.12400055, 58.12400055, 72.15100098, 72.15100098, 86.17790222, 100.2050018, 114.2320023, 18.01510048, 32.04190063]  # молярная масса, кг/кмоль (г/моль)
        Antoine = [[31.35, -1307.52, 0.0, -3.26134, 2.94e-05, 2.0, 0.0], [44.0103, -2568.82, 0.0, -4.97635, 1.46e-05, 2.0, 0.0], [52.3785, -3490.55, 0.0, -6.10875, 1.12e-05, 2.0, 0.0], [58.7845, -4136.68, 0.0, -7.01666, 1.04e-05, 2.0, 0.0], [66.945, -4604.09, 0.0, -8.25491, 1.16e-05, 2.0, 0.0], [66.7563, -5059.18, 0.0, -8.08935, 9.25e-06, 2.0, 0.0], [63.3315, -5117.78, 0.0, -7.48305, 7.77e-06, 2.0, 0.0], [70.4265, -6055.6, 0.0, -8.37865, 6.62e-06, 2.0, 0.0], [78.3285, -6947.0, 0.0, -9.44866, 6.47e-06, 2.0, 0.0], [86.997, -7890.6, 0.0, -10.6255, 6.47e-06, 2.0, 0.0], [65.9278, -7227.53, 0.0, -7.17695, 4.03e-06, 2.0, 0.0], [59.8373, -6282.89, 0.0, -6.37873, 4.62e-06, 2.0, 0.0]]  # набор коэффициентов Антуана для расчета полинома для каждого к-та
        CAS = ['74-82-8', '74-84-0', '74-98-6', '75-28-5', '106-97-8', '78-78-4', '109-66-0', '110-54-3', '142-82-5', '111-65-9', '7732-18-5', '67-56-1']  # импорт CAS номеров компонентов
        Hform = [-74900.0, -84738.0, -103890.0, -134590.0, -126190.0, -154590.0, -146490.0, -167290.0, -187890.0, -208590.0, -241814.0, -201290.0]  # Стандартная энтальпия образования вещества при 1 атм и 298.15 К, Дж/моль
        a_H = [[-12.98, 2.36459, -0.00213247, 5.66e-06, -3.72e-09, 8.61e-13, 1.0, 0.0], [-1.7675, 1.1429, -0.0003236, 4.24e-06, -3.39e-09, 8.82e-13, 1.0, 0.0], [39.4889, 0.395, 0.00211409, 3.96e-07, -6.67e-10, 1.68e-13, 1.0, 0.0], [30.903, 0.1533, 0.00263479, 7.27e-08, -7.28e-10, 2.37e-13, 1.0, 0.0], [67.721, 0.00854058, 0.00327699, -1.11e-06, 1.77e-10, -6.4e-15, 1.0, 0.0], [64.25, -0.131798, 0.003541, -1.33e-06, 2.51e-10, -1.3e-14, 1.0, 0.0], [63.198, -0.0117017, 0.0033164, -1.17e-06, 2e-10, -8.66e-15, 1.0, 0.0], [74.513, -0.096697, 0.00347649, -1.32e-06, 2.52e-10, -1.35e-14, 1.0, 0.0], [71.41, -0.0968949, 0.003473, -1.33e-06, 2.56e-10, -1.38e-14, 1.0, 0.0], [126.507, -0.2701, 0.00399829, -1.97e-06, 6.23e-10, -9.38e-14, 1.0, 0.0], [-5.7296, 1.9145, -0.00039574, 8.76e-07, -4.95e-10, 1.04e-13, 1.0, 0.0], [0.0, 0.6602, 0.0011072, 2.69e-07, -2.23e-10, 0.0, 1.0, 0.0]]  # набор коэффициентов энтальпии для расчета полинома для каждого к-та
        a_Cp = [[2.36459, 4.5474283009e-06, 1.81321496e-16, 1.9150131456e-34, 4.731684263603009e-61], [1.1429, 1.0471696e-07, 7.622502400000001e-17, 1.3206836241000002e-34, 5.3375619130243196e-61], [0.395, 4.469376528099999e-06, 6.2099136e-20, 1.9792622232099997e-37, 1.3382782156799996e-64], [0.1533, 6.9421183441e-06, 3.842405829999999e-22, 2.8088304025599996e-37, 7.477247049569999e-64], [0.00854058, 1.0738663460100001e-05, -1.3676309999999998e-18, 9.815062410000002e-40, -1.073741824e-71], [-0.131798, 1.2538680999999999e-05, -2.3526369999999998e-18, 3.969126001000001e-39, -3.71293e-70], [-0.0117017, 1.0998508960000002e-05, -1.601613e-18, 1.6000000000000002e-39, -4.870678456765758e-71], [-0.096697, 1.20859827201e-05, -2.2999680000000003e-18, 4.032758016000001e-39, -4.4840334374999985e-70], [-0.0968949, 1.2061729e-05, -2.3526369999999998e-18, 4.294967296e-39, -5.0049003168000005e-70], [-0.2701, 1.59863229241e-05, -7.645373000000002e-18, 1.5064412064100003e-37, -7.26129685547168e-66], [1.9145, 1.566101476e-07, 6.722213759999999e-19, 6.003725062500002e-38, 1.2166529024e-65], [0.6602, 1.22589184e-06, 1.9465108999999997e-20, 2.472973441e-39, 0.0]] # набор коэффициентов энтальпии (дифф-л) для расчета полинома для каждого к-та

        kij = [
            [0.0, 0.00224, 0.00683, 0.01311, 0.0123, 0.01763, 0.01793, 0.02347, 0.02886, 0.03416, 0.5, -0.035], 
            [0.00224, 0.0, 0.00126, 0.00457, 0.0041, 0.00741, 0.0078, 0.01141, 0.01532, 0.01932, 0.5, 0.045], 
            [0.00683, 0.00126, 0.0, 0.00104, 0.00082, 0.00258, 0.0027, 0.00514, 0.00789, 0.01085, 0.48, 0.06], 
            [0.01311, 0.00457, 0.00104, 0.0, 1e-05, 0.00035, 0.00039, 0.00157, 0.00322, 0.00521, 0.48, 0.069], 
            [0.0123, 0.0041, 0.00082, 1e-05, 0.0, 0.0005, 0.00055, 0.00187, 0.00365, 0.00575, 0.48, 0.069], 
            [0.01763, 0.00741, 0.00258, 0.00035, 0.0005, 0.0, 0.0, 0.00044, 0.00146, 0.00288, 0.5, 0.06], 
            [0.01793, 0.0078, 0.0027, 0.00039, 0.00055, 0.0, 0.0, 0.00039, 0.0074, 0.00276, 0.48, 0.06], 
            [0.02347, 0.01141, 0.00514, 0.00157, 0.00187, 0.00044, 0.00039, 0.0, -0.0078, 0.00107, 0.5, 0.051], 
            [0.02886, 0.01532, 0.00789, 0.00322, 0.00365, 0.00146, 0.0074, -0.0078, 0.0, 0.00024, 0.5, 0.08], 
            [0.03416, 0.01932, 0.01085, 0.00521, 0.00575, 0.00288, 0.00276, 0.00107, 0.00024, 0.0, 0.48, 0.09], 
            [0.5, 0.5, 0.48, 0.48, 0.48, 0.5, 0.48, 0.5, 0.5, 0.48, 0.0, -0.18], 
            [-0.035, 0.045, 0.06, 0.069, 0.069, 0.06, 0.06, 0.051, 0.08, 0.09, -0.18, 0.0]
            ]
        
        self.comp_props = {}
        for c, pc, tc, wc, mrc, antc, casc, hc, ahc, acpc, ki in zip(all_comps, p_cr, T_cr, w, Mr, Antoine, CAS, Hform, a_H, a_Cp, kij):
            self.comp_props.update({c : {
                                    'p_cr': pc,
                                    'T_cr': tc,
                                    'w': wc,
                                    'Mr': mrc,
                                    'Antoine': antc,
                                    'CAS': casc,
                                    'Hform': hc,
                                    'a_H': ahc,
                                    'a_Cp': acpc}})
            ckij = {}
            for cj, kj in zip(all_comps, ki):
                ckij.update({cj : kj})
            self.comp_props[c].update({'kij': ckij})
        
        self.p_cr = []  # критические давления, кПа
        self.T_cr = []  # критические температуры, С
        self.w = []  # ацентрический фактор
        self.Mr = []  # молярная масса, кг/кмоль (г/моль)
        self.Antoine = []  # набор коэффициентов Антуана для расчета полинома для каждого к-та
        self.CAS = []  # импорт CAS номеров компонентов
        self.Hform = []  # Стандартная энтальпия образования вещества при 1 атм и 298.15 К, Дж/моль
        self.a_H = []  # набор коэффициентов энтальпии для расчета полинома для каждого к-та
        self.a_Cp = []
        self.kij = []

        for comp in self.comps:
            self.p_cr.append(self.comp_props[comp]['p_cr'])
            self.T_cr.append(self.comp_props[comp]['T_cr'])
            self.w.append(self.comp_props[comp]['w'])
            self.Mr.append(self.comp_props[comp]['Mr'])
            self.Antoine.append(self.comp_props[comp]['Antoine'])
            self.CAS.append(self.comp_props[comp]['CAS'])
            self.Hform.append(self.comp_props[comp]['Hform'])
            self.a_H.append(self.comp_props[comp]['a_H'])
            self.a_Cp.append(self.comp_props[comp]['a_Cp'])
            ckij = []
            for cj in self.comps:
                ckij.append(self.comp_props[comp]['kij'][cj])
            self.kij.append(ckij)
        

class PengRobinson:

    R = 8.314
    # расходы газа, конденсата и воды, кг/ч

    def __init__(self, composition, pressure, temperature, mass_flow, conz_type='mol'):
        self.nm = list(composition.keys())
        #self.nm = ['Methane', 'Ethane', 'Propane', 'i-Butane', 'n-Butane', 'i-Pentane', 'n-Pentane', 'n-Hexane', 'n-Heptane', 'n-Octane', 'H2O', 'Methanol']
        self.N_comp = len(self.nm)
        self.p = pressure  # в кПа
        self.T = temperature  # в цельсиях
        self.G = mass_flow  # в кг/ч

        self.polar = ['-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '+', '+']

         # базово полярные (polar) всегда вода и спирты, как и любая ароматика, растворяющаяся в воде (напр. Фенол). Для
        # таких веществ ставится "+", чтобы константы считались по Антуану

        self.conz_type = conz_type  # mol или mass
        self.concec = []
        self.polar = []
        for c in composition:
            self.concec.append(composition[c])
            if c in ['H2O', 'Methanol']:
                self.polar.append('+')
            else:
                self.polar.append('-')
        #self.concec = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0]
        self.sep_type = '3'  # тип сепаратора: '2' - двухфазный, '3' - трехфазный

        self.mw_avg = 0

        if self.conz_type == 'mol':
            self.conz = self.concec
            self.w_conz = [0 for i in range(self.N_comp)]
        elif self.conz_type == 'mass':
            self.w_conz = self.concec
            self.conz = [0 for i in range(self.N_comp)]
        else:
            print("Неправильно указана размерность концентрации (mol или mass)")

        # мольные доли компонентов газа, конденсата и воды, кг/ч
        self.conz_y = [0 for i in range(self.N_comp)]
        self.conz_x = [0 for i in range(self.N_comp)]
        self.conz_w = [0 for i in range(self.N_comp)]

        # константы равновесия
        self.K = [0 for i in range(self.N_comp)]  # по Вильсону
        self.Kb = [0 for i in range(self.N_comp)]  #  по Антуану

        # константы компонентов для уравнения PR из базы SQL
        properties = ComponentProperties(self.nm)
        self.p_cr = properties.p_cr
        self.T_cr = properties.T_cr
        self.w = properties.w
        self.Mr = properties.Mr
        self.Antoine = properties.Antoine
        self.CAS = properties.CAS
        self.Hform = properties.Hform
        self.a_H = properties.a_H
        self.a_Cp = properties.a_Cp
        self.kij = properties.kij

        # величины для пересчета, фугитивности, коэффициенты уравнения PR
        self.p_cr_am = [0 for i in range(self.N_comp)]  # критические давления, Па
        self.T_cr_am = [0 for i in range(self.N_comp)]  # критические температуры, К
        self.p_r = [0 for i in range(self.N_comp)]  # приведенные давления
        self.T_r = [0 for i in range(self.N_comp)]  # приведенный температуры
        self.Ab = [[0 for i in range(self.N_comp)] for j in range(self.N_comp)]  # aij в уравнении PR
        self.Bp = [0 for i in range(self.N_comp)]  # bi в уравнении PR
        self.Ap = [0 for i in range(self.N_comp)]  # ai в уравнении PR
        self.Fiv = [0 for i in range(self.N_comp)]  # коэф-т фугитивности компонентов паровой фазы
        self.Fil = [0 for i in range(self.N_comp)]  # коэф-т фугитивности компонентов жидких фаз
        self.fiv = [0 for i in range(self.N_comp)]   # фугитивность компонентов паровой фазы
        self.fil = [0 for i in range(self.N_comp)]   # фугитивность компонентов жидких фаз
        self.OldK = [0 for i in range(self.N_comp)]  # старые константы равновесия по Вильсону
        self.OldKb = [0 for i in range(self.N_comp)]  # старые константы равновесия по Антуану
        self.kij = [[0 for i in range(self.N_comp)] for j in range(self.N_comp)]  # бинарные коэффициенты
        self.m = [0 for i in range(self.N_comp)]  # mi в уравнении PR
        self.alp = [0 for i in range(self.N_comp)]  # alfa(w, T) i в уравнении PR

        self.Gv = 0
        self.Gl = 0
        self.Gw = 0

        self.p_am = 0  # давление, Па
        self.T_am = 0  # температура, К

        self.Av = 0  # коэффициент a ур-я PR для газовой фазы
        self.Bv = 0  # коэффициент b ур-я PR для газовой фазы
        self.Al = 0  # -..- для жидких фаз
        self.Bl = 0  # -..- для жидких фаз
        self.e = 0  # доля отгона газа (vapour fraction)
        self.Zv = 0  # Коэффициент сжимаемости газовой фазы
        self.Zl = 0  # Коэффициент сжимаемости жидкой фазы
        self.Gm = 0  # Мольный расход потока, кмоль/ч


        self.GH2O = 0  # вода, которая удаляется на первом этапе, кмоль/ч
        self.step = ''  # нужно для понимания, по какой методике считать погрешность расчета - через константы или фугитивности

        # набор параметров уравнения PR для расчета энтальпии, сохраняются, так как рассчитываются на разных этапах
        self.H_y = []
        self.H_x = []
        self.H_w = []
        self.H = 0  # энтальпия потока
        self.S = 0  # энтропия потока
        self.Cp = 0 # массовая теплоемкость потока, Дж/кг*K
        self.Cv = 0  # массовая теплоемкость потока, Дж/кг*K
        self.y = 0  # показатель адиабаты
        self.v = 0 # молярный объем газа (потока), м3/моль
        self.ro = 0 # плотность при рабочих условиях, г/см3

        # словари для записей параметров потоков
        self.gas = {}
        self.liq = {}
        self.water = {}
        self.ishod_sm = {}  # состав исходного потока

        self.reid_pressure = None
        self.hc_dew_point = None
        self.water_dew_point = None

    def calculation(self):

        if self.conz_type == 'mass':
            self.calc_conz_w_to_x()
        elif self.conz_type == 'mol':
            pass

        self.mw_avg = self.get_Mr()

        self.ishod_sm['zi'] = [self.conz[i] for i in range(len(self.conz))]
        
        # self.reid_pressure = self.reid_pressure_calc()

        if 'H2O' in self.nm:
            if self.conz[self.nm.index('H2O')] != 1:
                self.H2O_delete()
            else:
                pass
        else:
            self.Gm = self.G / self.get_Mr()

        self.step = '1'
        self.pred_calculation()
        self.PR_calc()        

        self.proverka_phase()
        self.ro = self.ro_calc()
        self.reid_pressure = self.reid_pressure_calc()
        self.hc_dew_point = self.hc_dew_point_calc()
        self.water_dew_point = self.water_dew_point_calc()
        self.H = self.H_calc()
        # print(self.H)
        # self.S = self.S_calc()
        self.Cp = self.Cp_mass_calc()[0]
        self.Cv = self.Cp_mass_calc()[1]
        self.y = self.Cp_mass_calc()[2]

    def PR_calc(self):
        self._calculation()
        while self.calc_pogr() > 1e-2:
            self.memoryK()
            if not self._calculation():
                break

    def calc_conz_w_to_x(self):
        '''Пересчет из массовых долей в мольные'''
        w = 0
        for i in range(self.N_comp):
            w += self.w_conz[i] / self.Mr[i]
        for i in range(self.N_comp):
            self.conz[i] = (self.w_conz[i] / self.Mr[i]) / w

    def H2O_delete(self):
        '''Удаление воды и пересчет концентраций в потоке до 1'''
        self.Gm = self.G / self.get_Mr()
        Gm_with_water = self.Gm
        a = self.nm.index('H2O')
        self.GH2O = self.Gm * self.conz[a]

        Gm_new = self.Gm - self.Gm * self.conz[a]
        self.Gm = Gm_new
        self.conz[a] = 0
        for i in range(self.N_comp):
            self.conz[i] = (self.conz[i] * Gm_with_water) / Gm_new

    def get_Mr(self):
        '''Функция для расчета средней молекулярной массы'''
        Mr_s = 0
        for i in range(self.N_comp):
            Mr_s += self.conz[i] * self.Mr[i]
        return Mr_s

    def pred_calculation(self):
        self.recalc_P_T()
        self.calcK_Wilson()
        self.memoryK()

    def recalc_P_T(self):
        '''Пересчет температуры в Кельвины, давления в Паскали'''
        self.p_am = self.p * 1000
        self.T_am = self.T + 273.15
        for i in range(self.N_comp):
            self.p_cr_am[i] = self.p_cr[i] * 1000
            self.T_cr_am[i] = self.T_cr[i] + 273.15
            self.p_r[i] = self.p_am / self.p_cr_am[i]
            self.T_r[i] = self.T_am / self.T_cr_am[i]

    def calcK_Wilson(self):
        '''Расчет констант равновесия по уравнению Вильсона для начального приближения'''
        for i in range(self.N_comp):
            self.K[i] = (1 / self.p_r[i]) * math.exp(
                5.3727 * (1 + self.w[i]) * (1 - (1 / self.T_r[i])))

    def memoryK(self):
        for i in range(self.N_comp):
            self.OldK[i] = self.K[i]

    def _calculation(self):
        def func(Z1, Z2, Z3):
            if Z1 != 'null':
                return Z1
            if Z2 != 'null':
                return Z2
            if Z3 != 'null':
                return Z3

        self.calc_x_y()
        self.calc_Ab_Bp()
        self.calc_Av_Bv()
        # куб уравнение

        Av = self.Av * (self.p_am / (self.T_am ** 2 * self.R ** 2))  # параметр А для газа в кубическом уравнении PR
        Bv = self.Bv * (self.p_am / (self.R * self.T_am))  # параметр B для газа в кубическом уравнении PR

        self.calc_Al_Bl()
        # куб уравнение

        Al = self.Al * (self.p_am / (self.T_am ** 2 * self.R ** 2))  # параметр А для жидкостей в кубическом уравнении PR
        Bl = self.Bl * (self.p_am / (self.R * self.T_am))  # параметр В для жидкостей в кубическом уравнении PR

        Z1, Z2, Z3 = self.calcCubeExp(Bv - 1, Av - 2 * Bv - 3 * math.pow(Bv, 2),
                                      (-Av + math.pow(Bv, 2) + Bv) * Bv)

        # Zv должно быть максимумом из трех корней, поэтому все null я превращаю в 0
        self.Zv = func(Z1, Z2, Z3)
        if Z1 == 'null':
            Z1 = 0
        if Z2 == 'null':
            Z2 = 0
        if Z3 == 'null':
            Z3 = 0
        self.Zv = max(Z1, Z2, Z3)


        Z1, Z2, Z3 = self.calcCubeExp(Bl - 1, Al - 2 * Bl - 3 * math.pow(Bl, 2),
                                      (-Al + math.pow(Bl, 2) + Bl) * Bl)

        # Zl должно быть минимумом из трех корней, поэтому все null я превращаю в 1
        self.Zl = func(Z1, Z2, Z3)
        if Z1 == 'null':
            Z1 = 1
        if Z2 == 'null':
            Z2 = 1
        if Z3 == 'null':
            Z3 = 1
        self.Zl = min(Z1, Z2, Z3)

        # Z в обоих случаях должна быть больше 0 и меньше 1. 4 строки ниже для того, чтобы если минимумом из трех для
        # Zl окажется 0, расчет тормозился
        if self.Zl == 0:
            self.Zl = 1
        else:
            pass

        # Zv и Zl не должны быть равны 1. 1 говорит о том, что доля отгона 1 или 0, соответственно нет смысла считать
        # фугитивности, состав фаз уже определен
        if (self.Zv == 1):
            return False
        if (self.Zl == 1):
            return False


        self.calc_Fiv_Fil()
        self.calcK_PR()
        return True

    def calc_x_y(self):
        def FazeCalc():
            ''' расчет на фазовое состояние (одна или две фазы)0- две фазы, 1- 1 одна жидкая, 2- одна газовая, 3 - точка начала кипения, 4- точка начала росы'''
            S1 = 0
            S2 = 0
            for i in range(self.N_comp):
                S1 += self.conz[i] * self.K[i]
                if self.K[i] > 0:
                    S2 += self.conz[i] / self.K[i]
            # print(S1, S2)
            if S1 > 1 and S2 > 1:
                return 0  # 2 фазы
            if S1 < 1:
                return 1  # жидкая
            if S2 < 1:
                return 2  # газовая
            if S1 == 1:
                return 3  # точка кипения
            if S2 == 1:
                return 4  # точка росы
            return -1

        def Rashford_Rice(_fc):  # расчет мольной доли отгона
            def fn(_e):
                s = 0
                for i in range(self.N_comp):
                    if self.K[i] > 1 * 10 ** (-15):
                        s += self.conz[i] * (self.K[i] - 1) / (1 + _e * (self.K[i] - 1))
                return s

            eps = 1e-9
            a = 0
            b = 1
            _e = (a + b) / 2
            if _fc == 3:
                _e = 0
            if _fc == 4:
                _e = 1
            if fn(a) * fn(b) < 0:
                while abs(a - b) > eps:
                    _e = (a + b) / 2
                    if fn(a) * fn(_e) > 0:
                        a = _e
                    else:
                        b = _e
                    _e = (a + b) / 2

            for i in range(self.N_comp):
                self.conz_x[i] = self.conz[i] / (1 + _e * (self.K[i] - 1))
                self.conz_y[i] = self.K[i] * self.conz_x[i]

            self.e = _e

        fc = FazeCalc()
        e_old = self.e
        if fc in [0, 3, 4]:
            Rashford_Rice(fc)
            if self.e >= 1:
                self.e = 1
                for i in range(self.N_comp):
                    self.conz_x[i] = 0
                    self.conz_y[i] = self.conz[i]

            if self.e <= 0:
                self.e = 0
                for i in range(self.N_comp):
                    self.conz_x[i] = self.conz[i]
                    self.conz_y[i] = 0

        if fc == 1:
            self.e = 0
            for i in range(self.N_comp):
                self.conz_x[i] = self.conz[i]
                self.conz_y[i] = 0

        if fc == 2:
            self.e = 1
            for i in range(self.N_comp):
                self.conz_x[i] = 0
                self.conz_y[i] = self.conz[i]

    def calc_Ab_Bp(self):
        '''Расчет вириальных коэффициентов'''
        for i in range(self.N_comp):
            self.Bp[i] = (0.077796074 * self.R * self.T_cr_am[i] / self.p_cr_am[i])
            if self.w[i] < 0.491:
                self.m[i] = 0.37464 + 1.54226 * self.w[i] - 0.26992 * self.w[i] * self.w[i]
            else:
                self.m[i] = 0.379642 + 1.48503 * self.w[i] - 0.164423 * self.w[i] ** 2 + 0.016666 * self.w[i] ** 3
            self.alp[i] = math.pow(1 + self.m[i] * (1 - math.sqrt(self.T_r[i])), 2)
            self.Ap[i] = (0.457235529 * self.alp[i] * self.R**2 * self.T_cr_am[i]**2/self.p_cr_am[i])

        for i in range(self.N_comp):
            for j in range(self.N_comp):
                # if (i!=j):
                self.Ab[i][j] = (1 - self.kij[i][j]) * math.sqrt(self.Ap[i] * self.Ap[j])

    def calc_Av_Bv(self):
        '''Расчет вириальных коэффициентов'''
        self.Av = 0
        for i in range(self.N_comp):
            for j in range(self.N_comp):
                self.Av += self.conz_y[i] * self.conz_y[j] * self.Ab[i][j]
        self.Bv = 0
        for i in range(self.N_comp):
            self.Bv += self.conz_y[i] * self.Bp[i]

    def calc_Al_Bl(self):
        self.Al = 0
        for i in range(self.N_comp):
            for j in range(self.N_comp):
                self.Al += self.conz_x[i] * self.conz_x[j] * self.Ab[i][j]
        self.Bl = 0
        for i in range(self.N_comp):
            self.Bl += self.conz_x[i] * self.Bp[i]

    def calcCubeExp(self, A, B, C):
        '''Решение кубического уравнения тригонометрической формулой Виета
        x^3+Ax^2+Bx+C=0
        '''

        def sgn(x):
            if x > 0:
                return 1
            if x == 0:
                return 0
            if x < 0:
                return -1

        x1 = 'null'
        x2 = 'null'
        x3 = 'null'
        Q = (A * A - 3 * B) / 9
        R = (2 * A * A * A - 9 * A * B + 27 * C) / 54
        S = Q * Q * Q - R * R
        if S > 0:
            fi = 1 / 3 * math.acos(R / math.pow(Q, 3 / 2))
            x1 = -2 * math.sqrt(Q) * math.cos(fi) - A / 3
            if x1 < 0:
                x1 = 'null'
            x2 = -2 * math.sqrt(Q) * math.cos(fi + 2 / 3 * 3.14) - A / 3
            if x2 < 0:
                x2 = 'null'
            x3 = -2 * math.sqrt(Q) * math.cos(fi - 2 / 3 * 3.14) - A / 3
            if x3 < 0:
                x3 = 'null'
        if S < 0:
            if Q > 0:
                fi = 1 / 3 * math.acosh(abs(R) / math.pow(Q, 3 / 2))
                x1 = -2 * sgn(R) * math.sqrt(Q) * math.cosh(fi) - A / 3
                if x1 < 0:
                    x1 = 'null'
            if Q < 0:
                sss = abs(R) / math.pow(abs(Q), 3 / 2)
                fi = 1 / 3 * math.asinh(sss)
                x1 = -2 * sgn(R) * math.sqrt(abs(Q)) * math.sinh(fi) - A / 3
                if x1 < 0:
                    x1 = 'null'
            if Q == 0:
                x1 = -math.pow(C - A * A * A / 27, 1 / 3) - A / 3
                if x1 < 0:
                    x1 = 'null'
        if S == 0:
            x1 = -2 * sgn(R) * math.sqrt(Q) - A / 3
            if x1 < 0:
                x1 = 'null'
            x2 = sgn(R) * math.sqrt(Q) - A / 3
            if x2 < 0:
                x2 = 'null'
        return x1, x2, x3

    def calc_Fiv_Fil(self):
        for i in range(self.N_comp):
            sum_x = 0
            sum_y = 0
            for j in range(self.N_comp):
                sum_x += self.Ab[i][j] * self.conz_x[j]
                sum_y += self.Ab[i][j] * self.conz_y[j]

            Av = self.Av * (self.p_am / (self.T_am ** 2 * self.R ** 2))
            Bv = self.Bv * (self.p_am / (self.R * self.T_am))
            Al = self.Al * (self.p_am / (self.T_am ** 2 * self.R ** 2))
            Bl = self.Bl * (self.p_am / (self.R * self.T_am))


            self.Fiv[i] = math.exp((self.Zv - 1) * self.Bp[i] / self.Bv - math.log(self.Zv - Bv) - \
                                   (Av / (2.828 * Bv)) * ((2 * sum_y / self.Av) - (self.Bp[i] / self.Bv)) * math.log(
                (self.Zv + (1 + math.sqrt(2)) * Bv) / (self.Zv + (1 - math.sqrt(2)) * Bv)))

            self.Fil[i] = math.exp((self.Zl - 1) * self.Bp[i] / self.Bl - math.log(self.Zl - Bl) - \
                                   (Al / (2.828 * Bl)) * ((2 * sum_x / self.Al) - (self.Bp[i] / self.Bl)) * math.log(
                (self.Zl + (1 + math.sqrt(2)) * Bl) / (self.Zl + (1 - math.sqrt(2)) * Bl)))

            self.fil[i] = self.Fil[i]*self.conz_x[i]*self.p_am
            self.fiv[i] = self.Fiv[i]*self.conz_y[i]*self.p_am

    def calcK_PR(self):
        self.memoryK()
        for i in range(self.N_comp):
            self.K[i] = self.Fil[i] / self.Fiv[i]

    def calc_pogr(self):
        s = 0
        for i in range(self.N_comp):
            if self.step == '1':
                s += abs(self.OldK[i] - self.K[i])
            elif self.step == '2':
                if self.fiv[i] != 0:
                    s += abs(self.fil[i] / self.fiv[i] - 1)
                else:
                    s += 0
            else:
                print("Не задан этап расчета при расчете погрешности")
        return s

    def proverka_phase(self):
        '''Функция для проверки, сколько фаз существует в системе (выполняется после расчета первого двухфазного
        равновесия без воды)'''

        self.step = '2'

        self.H_y = [self.Av, self.Bv, self.Zv]
        self.H_x = [self.Al, self.Bl, self.Zl]


        self.gas['Gm, kmole/h'] = self.e * self.Gm
        self.liq['Gm, kmole/h'] = (1 - self.e) * self.Gm
        self.gas['t, C'] = self.T
        self.gas['p, kPa'] = self.p
        self.liq['t, C'] = self.T
        self.liq['p, kPa'] = self.p
        self.gas['n gas'] = self.gas['Gm, kmole/h'] / self.Gm
        self.liq['n liq'] = 1 - self.gas['n gas']
        self.water['Gm, kmole/h'] = 0
        self.water['t, C'] = self.T
        self.water['p, kPa'] = self.p
        self.water['n water'] = 0

        Mr_y = 0
        Mr_x = 0
        for i in range(self.N_comp):
            Mr_y += self.conz_y[i] * self.Mr[i]
            Mr_x += self.conz_x[i] * self.Mr[i]

        self.gas['Mr, kg/kmole'] = Mr_y
        self.liq['Mr, kg/kmole'] = Mr_x

        self.gas['G, kg/h'] = self.gas['Gm, kmole/h'] * Mr_y
        self.liq['G, kg/h'] = self.liq['Gm, kmole/h'] * Mr_x

        self.gas['ro, kg/m3'] = ((self.p * 1000) / (self.Zv * self.R * self.T_am)) / 1000 * Mr_y
        self.liq['ro, kg/m3'] = ((self.p * 1000) / (self.Zl * self.R * self.T_am)) / 1000 * Mr_x

        for i in range(self.N_comp):
            self.gas[self.nm[i]] = self.conz_y[i]
            self.liq[self.nm[i]] = self.conz_x[i]

        if self.GH2O == 0:
            pass
        elif self.GH2O != 0:
            if 0 < self.e < 1:
                self.calculation_gw()
                if self.e == 1:
                    pass
                elif 0 < self.e < 1:
                    self.H_y = [self.Av, self.Bv, self.Zv]
                    self.H_w = [self.Al, self.Bl, self.Zl]

                    self.gas['Gm, kmole/h'] = self.e * self.Gm
                    self.water['Gm, kmole/h'] = (1 - self.e) * self.Gm
                    self.gas['t, C'] = self.T
                    self.gas['p, kPa'] = self.p
                    self.water['t, C'] = self.T
                    self.water['p, kPa'] = self.p
                    self.gas['n gas'] = self.gas['Gm, kmole/h'] / self.Gm
                    self.water['n water'] = 1 - self.gas['n gas']


                    Mr_y = 0
                    Mr_w = 0
                    for i in range(self.N_comp):
                        Mr_y += self.conz_y[i] * self.Mr[i]
                        Mr_w += self.conz_x[i] * self.Mr[i]

                    self.gas['Mr, kg/kmole'] = Mr_y
                    self.water['Mr, kg/kmole'] = Mr_w

                    self.gas['G, kg/h'] = self.gas['Gm, kmole/h'] * Mr_y
                    self.water['G, kg/h'] = self.water['Gm, kmole/h'] * Mr_w

                    self.gas['ro, kg/m3'] = ((self.p * 1000) / (self.Zv * self.R * self.T_am)) / 1000 * Mr_y
                    self.water['ro, kg/m3'] = ((self.p * 1000) / (self.Zl * self.R * self.T_am)) / 1000 * Mr_w

                    for i in range(self.N_comp):
                        self.gas[self.nm[i]] = self.conz_y[i]
                        self.water[self.nm[i]] = self.conz_x[i]

                    self.calculation_nw()
                    if self.e == 1:
                        pass
                    elif 0 < self.e < 1:
                        self.H_x = [self.Av, self.Bv, self.Zv]
                        self.H_w = [self.Al, self.Bl, self.Zl]

                        self.liq['Gm, kmole/h'] = self.e * self.Gm
                        self.water['Gm, kmole/h'] = (1 - self.e) * self.Gm
                        self.liq['t, C'] = self.T
                        self.liq['p, kPa'] = self.p
                        self.water['t, C'] = self.T
                        self.water['p, kPa'] = self.p
                        self.liq['n liq'] = self.liq['Gm, kmole/h'] / self.Gm
                        self.water['n water'] = 1 - self.liq['n liq']

                        Mr_x = 0
                        Mr_w = 0
                        for i in range(self.N_comp):
                            Mr_x += self.conz_y[i] * self.Mr[i]
                            Mr_w += self.conz_x[i] * self.Mr[i]

                        self.liq['Mr, kg/kmole'] = Mr_x
                        self.water['Mr, kg/kmole'] = Mr_w

                        self.liq['G, kg/h'] = self.liq['Gm, kmole/h'] * Mr_x
                        self.water['G, kg/h'] = self.water['Gm, kmole/h'] * Mr_w

                        self.liq['ro, kg/m3'] = ((self.p * 1000) / (self.Zv * self.R * self.T_am)) / 1000 * Mr_x
                        self.water['ro, kg/m3'] = ((self.p * 1000) / (self.Zl * self.R * self.T_am)) / 1000 * Mr_w

                        for i in range(self.N_comp):
                            self.liq[self.nm[i]] = self.conz_y[i]
                            self.water[self.nm[i]] = self.conz_x[i]
                    elif self.e == 0:
                        self.liq.clear()
                        self.liq['n liq'] = 0
                        self.H_w = [self.Al, self.Bl, self.Zl]

                        self.water['Gm, kmole/h'] = (1 - self.e) * self.Gm
                        self.water['t, C'] = self.T
                        self.water['p, kPa'] = self.p
                        self.water['n water'] = 1

                        Mr_w = 0
                        for i in range(self.N_comp):
                            Mr_w += self.conz_x[i] * self.Mr[i]

                        self.water['Mr, kg/kmole'] = Mr_w
                        self.water['G, kg/h'] = self.water['Gm, kmole/h'] * Mr_w
                        self.water['ro, kg/m3'] = ((self.p * 1000) / (self.Zl * self.R * self.T_am)) / 1000 * Mr_w

                        for i in range(self.N_comp):
                            self.water[self.nm[i]] = self.conz_x[i]

            if self.e == 1:
                self.calculation_gw()
                if self.e == 1:
                    self.H_y = [self.Av, self.Bv, self.Zv]
                    self.gas['Gm, kmole/h'] = self.e * self.Gm
                    self.gas['t, C'] = self.T
                    self.gas['p, kPa'] = self.p
                    self.gas['n gas'] = self.gas['Gm, kmole/h'] / self.Gm

                    Mr_y = 0
                    for i in range(self.N_comp):
                        Mr_y += self.conz_y[i] * self.Mr[i]


                    self.gas['Mr, kg/kmole'] = Mr_y


                    self.gas['G, kg/h'] = self.gas['Gm, kmole/h'] * Mr_y


                    self.gas['ro, kg/m3'] = ((self.p * 1000) / (self.Zv * self.R * self.T_am)) / 1000 * Mr_y


                    for i in range(self.N_comp):
                        self.gas[self.nm[i]] = self.conz_y[i]
                    pass
                elif 0 < self.e < 1:
                    self.H_y = [self.Av, self.Bv, self.Zv]
                    self.H_w = [self.Al, self.Bl, self.Zl]

                    self.gas['Gm, kmole/h'] = self.e * self.Gm
                    self.water['Gm, kmole/h'] = (1 - self.e) * self.Gm
                    self.gas['t, C'] = self.T
                    self.gas['p, kPa'] = self.p
                    self.water['t, C'] = self.T
                    self.water['p, kPa'] = self.p
                    self.gas['n gas'] = self.gas['Gm, kmole/h'] / self.Gm
                    self.water['n water'] = 1 - self.gas['n gas']

                    Mr_y = 0
                    Mr_w = 0
                    for i in range(self.N_comp):
                        Mr_y += self.conz_y[i] * self.Mr[i]
                        Mr_w += self.conz_x[i] * self.Mr[i]

                    self.gas['Mr, kg/kmole'] = Mr_y
                    self.water['Mr, kg/kmole'] = Mr_w

                    self.gas['G, kg/h'] = self.gas['Gm, kmole/h'] * Mr_y
                    self.water['G, kg/h'] = self.water['Gm, kmole/h'] * Mr_w

                    self.gas['ro, kg/m3'] = ((self.p * 1000) / (self.Zv * self.R * self.T_am)) / 1000 * Mr_y
                    self.water['ro, kg/m3'] = ((self.p * 1000) / (self.Zl * self.R * self.T_am)) / 1000 * Mr_w

                    for i in range(self.N_comp):
                        self.gas[self.nm[i]] = self.conz_y[i]
                        self.water[self.nm[i]] = self.conz_x[i]

            if self.e == 0:
                self.calculation_nw()
                if self.e == 1:
                    self.H_x = [self.Av, self.Bv, self.Zv]

                    self.liq['Gm, kmole/h'] = self.e * self.Gm
                    self.liq['t, C'] = self.T
                    self.liq['p, kPa'] = self.p
                    self.liq['n liq'] = self.liq['Gm, kmole/h'] / self.Gm

                    Mr_x = 0
                    for i in range(self.N_comp):
                        Mr_x += self.conz_y[i] * self.Mr[i]

                    self.liq['Mr, kg/kmole'] = Mr_x

                    self.liq['G, kg/h'] = self.liq['Gm, kmole/h'] * Mr_x


                    self.liq['ro, kg/m3'] = ((self.p * 1000) / (self.Zv * self.R * self.T_am)) / 1000 * Mr_x

                    for i in range(self.N_comp):
                        self.liq[self.nm[i]] = self.conz_y[i]
                    pass
                elif 0 < self.e < 1:
                    self.H_x = [self.Av, self.Bv, self.Zv]
                    self.H_w = [self.Al, self.Bl, self.Zl]

                    self.liq['Gm, kmole/h'] = self.e * self.Gm
                    self.water['Gm, kmole/h'] = (1 - self.e) * self.Gm
                    self.liq['t, C'] = self.T
                    self.liq['p, kPa'] = self.p
                    self.water['t, C'] = self.T
                    self.water['p, kPa'] = self.p
                    self.liq['n liq'] = self.liq['Gm, kmole/h'] / self.Gm


                    Mr_x = 0
                    Mr_w = 0
                    for i in range(self.N_comp):
                        Mr_x += self.conz_y[i] * self.Mr[i]
                        Mr_w += self.conz_x[i] * self.Mr[i]

                    self.liq['Mr, kg/kmole'] = Mr_x
                    self.water['Mr, kg/kmole'] = Mr_w

                    self.liq['G, kg/h'] = self.liq['Gm, kmole/h'] * Mr_x
                    self.water['G, kg/h'] = self.water['Gm, kmole/h'] * Mr_w

                    self.liq['ro, kg/m3'] = ((self.p * 1000) / (self.Zv * self.R * self.T_am)) / 1000 * Mr_x
                    self.water['ro, kg/m3'] = ((self.p * 1000) / (self.Zl * self.R * self.T_am)) / 1000 * Mr_w

                    for i in range(self.N_comp):
                        self.liq[self.nm[i]] = self.conz_y[i]
                        self.water[self.nm[i]] = self.conz_x[i]
                elif self.e == 0:
                    self.liq.clear()
                    self.liq['n liq'] = 0
                    self.H_w = [self.Al, self.Bl, self.Zl]

                    self.water['Gm, kmole/h'] = (1 - self.e) * self.Gm
                    self.water['t, C'] = self.T
                    self.water['p, kPa'] = self.p
                    self.water['n water'] = 1

                    Mr_w = 0
                    for i in range(self.N_comp):
                        Mr_w += self.conz_x[i] * self.Mr[i]

                    self.water['Mr, kg/kmole'] = Mr_w
                    self.water['G, kg/h'] = self.water['Gm, kmole/h'] * Mr_w
                    self.water['ro, kg/m3'] = ((self.p * 1000) / (self.Zl * self.R * self.T_am)) / 1000 * Mr_w

                    for i in range(self.N_comp):
                        self.water[self.nm[i]] = self.conz_x[i]

        Ggas = self.gas['Gm, kmole/h'] if 'Gm, kmole/h' in self.gas else 0
        Gliq = self.liq['Gm, kmole/h'] if 'Gm, kmole/h' in self.liq else 0
        Gwater = self.water['Gm, kmole/h'] if 'Gm, kmole/h' in self.water else 0

        self.gas['n gas'] = Ggas/(Ggas + Gliq + Gwater)
        self.liq['n liq'] = Gliq / (Ggas + Gliq + Gwater)
        self.water['n water'] = Gwater / (Ggas + Gliq + Gwater)

        self.gas['G, kg/h'] = 0 if self.gas['n gas'] == 0 else self.gas['G, kg/h']
        self.liq['G, kg/h'] = 0 if self.liq['n liq'] == 0 else self.liq['G, kg/h']
        self.water['G, kg/h'] = 0 if self.water['n water'] == 0 else self.water['G, kg/h']

        if self.sep_type == '2' and self.liq['Gm, kmole/h'] != 0 and self.water['Gm, kmole/h']:
            self.liq['Mr, kg/kmole'] = self.liq['Mr, kg/kmole'] * self.liq['Gm, kmole/h'] / (self.liq['Gm, kmole/h'] + self.water['Gm, kmole/h']) + self.water['Mr, kg/kmole'] * self.water['Gm, kmole/h'] / (self.liq['Gm, kmole/h'] + self.water['Gm, kmole/h'])
            for i in range(self.N_comp):
                self.liq[self.nm[i]] = (self.liq[self.nm[i]] * self.liq['Gm, kmole/h'] + self.water[self.nm[i]] * self.water['Gm, kmole/h']) / (self.liq['Gm, kmole/h'] + self.water['Gm, kmole/h'])
            self.liq['ro, kg/m3'] = self.liq['ro, kg/m3'] * self.liq['G, kg/h'] / (self.liq['G, kg/h'] + self.water['G, kg/h']) + self.water['ro, kg/m3'] * self.water['G, kg/h'] / (self.liq['G, kg/h'] + self.water['G, kg/h'])
            self.liq['Gm, kmole/h'] = self.liq['Gm, kmole/h'] + self.water['Gm, kmole/h']
            self.liq['G, kg/h'] = self.liq['G, kg/h'] + self.water['G, kg/h']
            self.liq['n liq'] = self.liq['n liq'] + self.water['n water']
            self.water['Gm, kmole/h'] = 0

    def calculation_gw(self):
        "расчет газ+вода"
        Gm_new = self.gas['Gm, kmole/h'] + self.GH2O
        self.Gm = Gm_new
        for i in range(self.N_comp):
            self.conz[i] = (self.conz_y[i]*self.gas['Gm, kmole/h'])/Gm_new
        self.conz[self.nm.index('H2O')] = self.GH2O/Gm_new

        for i in range(self.N_comp):
            self.K[i] = 10 ** 6 * (self.p_r[i] / self.T_r[i])
        a = self.polar.count('+')
        n = [i for i in range(0, len(self.polar)) if self.polar[i] == '+']
        for i in range(a):
            T_am = self.T_am
            if n[i] == self.nm.index('H2O'):
                if self.T_am < 281:
                    T_am = 281
                elif self.T_am > 647.3:
                    T_am = 647.3
                else:
                    T_am = self.T_am
            self.K[n[i]] = (self.Antoine[n[i]][0] + self.Antoine[n[i]][1] / (T_am + self.Antoine[n[i]][2]) +
                            self.Antoine[n[i]][3] * math.log(T_am) + self.Antoine[n[i]][4] * T_am **
                            self.Antoine[n[i]][5]) / self.p_am

        self.memoryK()
        self.PR_calc()

    def calculation_nw(self):
        "расчет конденсат+вода"
        Gm_new = self.liq['Gm, kmole/h'] + self.water['Gm, kmole/h']
        self.Gm = Gm_new
        for i in range(self.N_comp):
            self.conz[i] = (self.liq[self.nm[i]] * self.liq['Gm, kmole/h'] + self.water[self.nm[i]] * self.water['Gm, kmole/h']) / Gm_new

        a = self.polar.count('+')
        n = [i for i in range(0, len(self.polar)) if self.polar[i] == '+']


        for i in range(self.N_comp):
            self.K[i] = 10 ** 6 * (self.p_r[i] / self.T_r[i])
        for i in range(a):
            T_am = self.T_am
            if n[i] == self.nm.index('H2O'):
                if self.T_am < 281:
                        T_am = 281
                elif self.T_am > 647.3:
                        T_am = 647.3
                else:
                        T_am = self.T_am
            self.K[n[i]] = (self.Antoine[n[i]][0] + self.Antoine[n[i]][1]/(T_am+self.Antoine[n[i]][2]) + self.Antoine[n[i]][3] * math.log(T_am) + self.Antoine[n[i]][4] * T_am**self.Antoine[n[i]][5])/self.p_am
        self.memoryK()
        self.PR_calc()

    def H_calc(self):
        H_H = []
        H_form = []
        H = []
        Hr = 0
        Mr = 0

        for i in range(self.N_comp):
            Mr = 0
            self.conz = [0 for i in range(self.N_comp)]
            self.conz[i] = 1.0
            self.conz_type = 'mol'
            self.PR_calc()
            if self.e == 1:

                Av_y = self.H_y[0]
                Bv_y = self.H_y[1]
                Zv_y = self.H_y[2]
                Av = Av_y * (self.p_am / (self.T_am ** 2 * self.R ** 2))
                Bv = Bv_y * (self.p_am / (self.R * self.T_am))

                Mr = self.Mr[i]

                H_H_a = self.R * self.T_cr_am[i] * (
                        self.T_r[i] * (Zv_y - 1) - 2.078 * (1 + self.m[i]) * math.sqrt(self.alp[i]) * math.log(
                    (Zv_y + 2.414 * Bv) / (Zv_y - 0.414 * Bv)))

                H_form_a = self.Hform[i]

                H_a = (self.a_H[i][0] + self.a_H[i][1] * self.T_am + self.a_H[i][2] * (
                        self.T + 273.15) ** 2 +
                       self.a_H[i][3] * self.T_am ** 3 + self.a_H[i][4] * self.T_am ** 4
                       + self.a_H[i][5] * self.T_am ** 5) - (
                              self.a_H[i][0] + self.a_H[i][1] * 298.15 + self.a_H[i][
                          2] * 298.15 ** 2 +
                              self.a_H[i][3] * 298.15 ** 3 + self.a_H[i][4] * 298.15 ** 4
                              + self.a_H[i][5] * 298.15 ** 5)
                H.append(H_a + H_form_a / Mr + H_H_a / Mr)
            elif self.e == 0:
                Av_x = self.H_x[0]
                Bv_x = self.H_x[1]
                Zv_x = self.H_x[2]
                Ax = Av_x * (self.p_am / (self.T_am ** 2 * self.R ** 2))
                Bx = Bv_x * (self.p_am / (self.R * self.T_am))

                Mr = self.Mr[i]

                H_H_a = self.R * self.T_cr_am[i] * (
                        self.T_r[i] * (Zv_x - 1) - 2.078 * (1 + self.m[i]) * math.sqrt(self.alp[i]) * math.log(
                    (Zv_x + 2.414 * Bx) / (Zv_x - 0.414 * Bx)))

                H_form_a = self.Hform[i]

                H_a = (self.a_H[i][0] + self.a_H[i][1] * self.T_am + self.a_H[i][2] * (
                        self.T + 273.15) ** 2 +
                       self.a_H[i][3] * self.T_am ** 3 + self.a_H[i][4] * self.T_am ** 4
                       + self.a_H[i][5] * self.T_am ** 5) - (
                              self.a_H[i][0] + self.a_H[i][1] * 298.15 + self.a_H[i][
                          2] * 298.15 ** 2 +
                              self.a_H[i][3] * 298.15 ** 3 + self.a_H[i][4] * 298.15 ** 4
                              + self.a_H[i][5] * 298.15 ** 5)
                H.append(H_a + H_form_a / Mr + H_H_a / Mr)
            else:
                print("kipit")

        self.conz = self.ishod_sm['zi']
        self.PR_calc()
        H_y = 0
        H_x = 0
        H_w = 0
        conz_y = []
        conz_x = []
        conz_w = []

        if self.gas['n gas'] != 0:
            for i in range(self.N_comp):
                conz_y.append(self.gas[self.nm[i]])
            for i in range(self.N_comp):
                H_y += H[i] * conz_y[i] * self.Mr[i] / self.gas['Mr, kg/kmole']
            ngasmass = self.gas['G, kg/h']/(self.liq['G, kg/h']+self.gas['G, kg/h']+self.water['G, kg/h'])
        else:
            ngasmass = 0

        if self.liq['n liq'] != 0:
            for i in range(self.N_comp):
                conz_x.append(self.liq[self.nm[i]])
            for i in range(self.N_comp):
                H_x += H[i] * conz_x[i] * self.Mr[i] / self.liq['Mr, kg/kmole']
            nliqmass = self.liq['G, kg/h'] / (self.liq['G, kg/h'] + self.gas['G, kg/h'] + self.water['G, kg/h'])
        else:
            nliqmass = 0

        if self.water['n water'] != 0:
            for i in range(self.N_comp):
                conz_w.append(self.water[self.nm[i]])
            for i in range(self.N_comp):
                H_w += H[i] * conz_w[i] * self.Mr[i] / self.water['Mr, kg/kmole']
            nwatermass = self.water['G, kg/h'] / (self.liq['G, kg/h'] + self.gas['G, kg/h'] + self.water['G, kg/h'])
        else:
            nwatermass = 0

        Hr = H_y * ngasmass + H_x * nliqmass + H_w * nwatermass
        #print("Энтальпия:", Hr)
        return Hr

    def S_calc(self):
        S_S = []
        S_form = []
        S = []
        Sr = 0
        Mr = 0

        for i in range(self.N_comp):
            Mr = 0
            self.conz = [0 for i in range(self.N_comp)]
            self.conz[i] = 1.0
            self.conz_type = 'mol'
            self.PR_calc()
            if self.e == 1:
                a = self.H_y[0]
                b = self.H_y[1]
                Z = self.H_y[2]
            elif self.e == 0:
                a = self.H_x[0]
                b = self.H_x[1]
                Z = self.H_x[2]

            A = a * (self.p_am / (self.T_am ** 2 * self.R ** 2))
            B = b * (self.p_am / (self.R * self.T_am))


            Mr = self.Mr[i]
            dadt = self.dadT_p(i)
            v = self.v_calc()
            S_S_a = self.R * (Z - B) - math.log(self.p_am/101.325) - A/(2**1.5 * b * self.R * self.T_am) * self.T_am/a * dadt * math.log((v+(2**0.5 + 1)*b)/(v+(2**0.5 - 1)*b))

            S_form_a = self.Hform[i]

            S_a = (self.a_Cp[i][0] + self.a_Cp[i][1] * self.T_am + self.a_Cp[i][2] * (
                        self.T + 273.15) ** 2 +
                       self.a_Cp[i][3] * self.T_am ** 3 + self.a_Cp[i][4] * self.T_am ** 4) - (
                              self.a_Cp[i][0] + self.a_Cp[i][1] * 298.15 + self.a_Cp[i][
                          2] * 298.15 ** 2 +
                              self.a_Cp[i][3] * 298.15 ** 3 + self.a_Cp[i][4] * 298.15 ** 4)
            S.append(S_S_a+S_a)


        # print("entropy", S)
        self.conz = self.ishod_sm['zi']
        self.PR_calc()
        S_y = 0
        S_x = 0
        S_w = 0
        conz_y = []
        conz_x = []
        conz_w = []

        if self.gas['n gas'] != 0:
            for i in range(self.N_comp):
                conz_y.append(self.gas[self.nm[i]])
            for i in range(self.N_comp):
                S_y += S[i] * conz_y[i] * self.Mr[i] / self.gas['Mr, kg/kmole']
            ngasmass = self.gas['G, kg/h']/(self.liq['G, kg/h']+self.gas['G, kg/h']+self.water['G, kg/h'])
        else:
            ngasmass = 0

        if self.liq['n liq'] != 0:
            for i in range(self.N_comp):
                conz_x.append(self.liq[self.nm[i]])
            for i in range(self.N_comp):
                S_x += S[i] * conz_x[i] * self.Mr[i] / self.liq['Mr, kg/kmole']
            nliqmass = self.liq['G, kg/h'] / (self.liq['G, kg/h'] + self.gas['G, kg/h'] + self.water['G, kg/h'])
        else:
            nliqmass = 0

        if self.water['n water'] != 0:
            for i in range(self.N_comp):
                conz_w.append(self.water[self.nm[i]])
            for i in range(self.N_comp):
                S_w += S[i] * conz_w[i] * self.Mr[i] / self.water['Mr, kg/kmole']
            nwatermass = self.water['G, kg/h'] / (self.liq['G, kg/h'] + self.gas['G, kg/h'] + self.water['G, kg/h'])
        else:
            nwatermass = 0

        Sr = S_y * ngasmass + S_x * nliqmass + S_w * nwatermass
        #print("Энтропия:", Sr)
        return Sr

    def Cp_mass_calc(self):
        Cp = []
        Cv = []



        for i in range(self.N_comp):
            self.p_am = self.p/100
            self.p_cr_am[i] = self.p_cr[i] / 100
            self.R = 83.14
            self.p_r[i] = self.p_am/self.p_cr_am[i]
            self.conz = [0 for i in range(self.N_comp)]
            self.conz[i] = 1.0
            self.conz_type = 'mol'
            self.PR_calc()

            if self.e == 1:
                B = self.Bv
                A = self.Av
                Z = self.H_y[2]
            else:
                B = self.Bl
                A = self.Al
                Z = self.H_x[2]
            
            Av = 0.45723553 * self.alp[i] * (self.p_r[i]/self.T_r[i]**2)  # параметр А для газа в кубическом уравнении PR
            Bv = 0.077796074 * (self.p_r[i]/self.T_r[i])  # параметр B для газа в кубическом уравнении PR

        
            Cp_ig = (self.a_H[i][1] + 2 * self.a_H[i][2] * (
                    self.T + 273.15) ** 1 +
                   3 * self.a_H[i][3] * self.T_am ** 2 + 4 * self.a_H[i][4] * self.T_am ** 3
                   + 5 * self.a_H[i][5] * self.T_am ** 4) * self.Mr[i]
            Cv_ig = Cp_ig - self.R/10
            Cv_R = ((self.T_am * self.d2adT2_v(i) / (B * 8**0.5)) * math.log((Z + Bv * (1 + 2**0.5)) / (Z + Bv * (1 - 2**0.5)))) / 10
            Cp_R = Cv_R + self.T_am * self.dpdT_V(i) * self.dvdT_p(i)/10 - self.R/10
            Cpg = Cp_ig + Cp_R
            Cvg = Cv_ig + Cv_R
            Cp.append(Cpg/self.Mr[i])
            Cv.append(Cvg/self.Mr[i])

        self.conz = self.ishod_sm['zi']
        self.R = 8.314
        self.recalc_P_T()
        self.PR_calc()
        Cp_y = 0
        Cp_x = 0
        Cp_w = 0
        Cv_y = 0
        Cv_x = 0
        Cv_w = 0
        conz_y = []
        conz_x = []
        conz_w = []

        if self.gas['n gas'] != 0:
            for i in range(self.N_comp):
                conz_y.append(self.gas[self.nm[i]])
            for i in range(self.N_comp):
                Cp_y += Cp[i] * conz_y[i] * self.Mr[i] / self.gas['Mr, kg/kmole']
                Cv_y += Cv[i] * conz_y[i] * self.Mr[i] / self.gas['Mr, kg/kmole']
            ngasmass = self.gas['G, kg/h'] / (self.liq['G, kg/h'] + self.gas['G, kg/h'] + self.water['G, kg/h'])
        else:
            ngasmass = 0

        if self.liq['n liq'] != 0:
            for i in range(self.N_comp):
                conz_x.append(self.liq[self.nm[i]])
            for i in range(self.N_comp):
                Cp_x += Cp[i] * conz_x[i] * self.Mr[i] / self.liq['Mr, kg/kmole']
                Cv_x += Cv[i] * conz_y[i] * self.Mr[i] / self.gas['Mr, kg/kmole']
            nliqmass = self.liq['G, kg/h'] / (self.liq['G, kg/h'] + self.gas['G, kg/h'] + self.water['G, kg/h'])
        else:
            nliqmass = 0

        if self.water['n water'] != 0:
            for i in range(self.N_comp):
                conz_w.append(self.water[self.nm[i]])
            for i in range(self.N_comp):
                Cp_w += Cp[i] * conz_w[i] * self.Mr[i] / self.water['Mr, kg/kmole']
                Cv_w += Cv[i] * conz_y[i] * self.Mr[i] / self.gas['Mr, kg/kmole']
            nwatermass = self.water['G, kg/h'] / (self.liq['G, kg/h'] + self.gas['G, kg/h'] + self.water['G, kg/h'])
        else:
            nwatermass = 0

        Cp_sm = Cp_y * ngasmass + Cp_x * nliqmass + Cp_w * nwatermass
        Cv_sm = Cv_y * ngasmass + Cv_x * nliqmass + Cv_w * nwatermass
        y = Cp_sm/Cv_sm
        return Cp_sm, Cv_sm, y

    def reid_pressure_calc(self):
        temp = 37.7 + 273

        vol_fracs = []
        for nc, xc, mc in zip(self.nm, self.ishod_sm['zi'], self.Mr):
            mass_fraction = xc * mc / self.mw_avg

            vf = PengRobinson({nc: 1}, self.p, self.T, mass_fraction*self.G)
            vf.step = '1'
            vf.pred_calculation()
            vf.PR_calc()        
            if vf.e == 1:
                z_vf = vf.Zv
            else:
                z_vf = vf.Zl
            density_vf = mc * self.p / (z_vf * 8.314 * self.T_am)
            vol_fracs.append(mass_fraction * density_vf / self.ro)
            del vf

        reid_indices = []
        for i, vf in zip(list(range(len(self.nm))), vol_fracs):
            ant = math.pow((self.Antoine[i][0] + self.Antoine[i][1] / (temp + self.Antoine[i][2]) + self.Antoine[i][3] * math.log(temp) + self.Antoine[i][4] * temp **self.Antoine[i][5]) / (6.8947 * 1000), 1.25)
            reid_indices.append(vf * math.pow((self.Antoine[i][0] + self.Antoine[i][1] / (temp + self.Antoine[i][2]) + self.Antoine[i][3] * math.log(temp) + self.Antoine[i][4] * temp **self.Antoine[i][5]) / (6.8947 * 1000), 1.25))

        return math.pow(sum(reid_indices), 1/1.25) * (6.8947 * 1000)

    def water_dew_point_calc(self):

        if self.gas['G, kg/h'] == 0:
            return None

        wp_comp = {}
        for c, y in zip(self.nm, self.conz):
            wp_comp.update({c: y})
        dpw = PengRobinson(wp_comp, self.p, -50, self.gas['G, kg/h'])
        dpw.step = '1'
        dpw.pred_calculation()
        dpw.PR_calc()
        while  dpw.e != 1:
            dpw.step = '1'
            dpw.pred_calculation()
            dpw.PR_calc()            
            dpw.T += 1
            # print(dpw.T)
        
        return dpw.T
    
    def hc_dew_point_calc(self):

        if self.gas['G, kg/h'] == 0:
            return None
        
        gas_comp = {}
        for c, y in zip(self.nm, self.conz_y):
            gas_comp.update({c: y})
        dp = PengRobinson(gas_comp, self.p, -200, self.gas['G, kg/h'])
        if 'H2O' in self.nm:
            if dp.conz[dp.nm.index('H2O')] != 1:
                dp.H2O_delete()
            else:
                pass
        else:
            dp.Gm = dp.G / dp.get_Mr()
        dp.step = '1'
        dp.pred_calculation()
        dp.PR_calc()
        while dp.e != 1:
            dp.step = '1'
            dp.pred_calculation()
            dp.PR_calc()
            dp.T += 1
        
        return dp.T
    
    '''
    def phase_envelope(self):
        def envelope_calc():
            self.step = '1'
            self.pred_calculation()
            self.PR_calc()
        def polar_delete():
            Gm = self.G / self.get_Mr()
            self.Gm = Gm
            a = self.polar.count('+')
            n = [i for i in range(0, len(self.polar)) if self.polar[i] == '+']
            for i in range(a):
                if n[i] == self.nm.index('H2O') or n[i] == self.nm.index('Methanol'):
                    self.Gm -= self.Gm * self.conz[n[i]]
                    self.conz[n[i]] = 0
            for i in range(self.N_comp):
                self.conz[i] = (self.conz[i] * Gm) / self.Gm

        self.p = 201
        self.T = -150
        p = []
        p1 = []
        T = []
        T1 = []
        ppr = []
        Tpr = []
        i = 0
        polar_delete()

        # расчет фазовой кривой
        while self.p < 12000:
            # envelope_calc()
            while self.T < 50:
                envelope_calc()
                self.T += 1
                # envelope_calc()
                if self.e != 0 and self.e != 1:
                    ppr.append(self.p)
                    Tpr.append(self.T)
                else:
                    pass

            # помещает при каждом давлении найденную минимальную и максимальную температуру в начало и конец списка
            # для правильного построения графика
            if len(ppr) != 0 and len(Tpr) != 0:
                p.insert(i, ppr[1])
                p.insert(len(p) - i, ppr[1])
                T.insert(i, min(Tpr))
                T.insert(len(T) - i, max(Tpr))
                i += 1
            else:
                pass

            if self.p < 6000:
                self.p += 500
                self.T = -150
            else:
                self.p += 100
                self.T = -150
            ppr.clear()
            Tpr.clear()


        # расчет критической точки
        self.T = -150
        self.p = 500
        Tcr = []
        pcr = []
        while self.T < 101:
            if len(pcr) != 0:
                break
            envelope_calc()
            while self.p < 20000:
                self.p += 200
                envelope_calc()
                if self.e == 1:
                    pcr.append(self.p)
                    Tcr.append(self.T)
                else:
                    pass

            self.p = 500
            self.T += 1

        p1.append(max(pcr))
        T1.append(max(Tcr))

        self.T = T1[0]
        self.p = p1[0]
        pt = p[T1.index(self.T)]
        self.p = pt
        T1.clear()
        p1.clear()
        T1.append(self.T)
        p1.append(self.p)

        plt.plot(T1, p1, 'bo', label='критическая точка')
        plt.plot(T, p, 'r', label='т. кипения')
        print(T)
        print(p)
        plt.grid()
        plt.legend()
        plt.xlabel('Температура, С')
        plt.ylabel('Давление, кПа')
        plt.show()
        '''

    def calc_drossel(self, p2):
        self.conz = self.ishod_sm['zi']
        H_old = self.H_calc()
        pdr = p2
        self.p = pdr
        self.calculation()
        T2 = self.T
        H = self.H_calc()
        while H > H_old:
            T2 -= 0.1
            self.T = T2
            self.calculation()
            H = self.H_calc()
        print(T2)
        return T2, pdr

    def Methanol_STO(self):
        conz1 = []  # Состав газа до точки защиты
        for i in range(self.N_comp):
            conz1.append(self.gas[self.nm[i]])
        Z1 = self.H_y[2]  # К-т сжимаемости газа до
        T1 = self.T  # Температура газа до точки защиты, С
        p1 = self.p  # Давление газа до точки защиты (дросселя, тдт, то), кПа

        Ggasstd = self.gas['Gm, kmole/h'] * 1000 * 24.465 / 1000  # ст. м3/ч при 25 С и 101.3 кПа
        if self.liq['n liq'] != 0:
            Mrliq1 = self.liq['Mr, kg/kmole']
            Gk1 = self.liq['G, kg/h'] / (Ggasstd / 1000)  # количество конденсата до точки защиты, кг/1000 м3
        else:
            Mrliq1 = 0
            Gk1 = 0

        if 'Methanol' in self.nm and conz1[self.nm.index('Methanol')] != 0:
            if self.water['n water'] != 0:
                water = []
                for i in range(self.N_comp):
                    water.append(self.water[self.nm[i]])
                X1 = ((water[self.nm.index('Methanol')] * self.water['Gm, kmole/h'] * 32) / (water[self.nm.index('Methanol')] * self.water['Gm, kmole/h'] * 32 + water[self.nm.index('H2O')] * self.water['Gm, kmole/h'] * 18))*100  # Концентрация метанола в ВМР, приходящего на узел, %mass
                Gmeth1 = self.water['G, kg/h']  # Расход ВМР, приходящего на узел, кг/ч
            else:
                X1 = 0  # Концентрация метанола в ВМР, приходящего на узел, %
                Gmeth1 = 0 # Расход ВМР, приходящего на узел, кг/ч


        comp_plotn_std = {'Methane': 0.55, 'Ethane': 0.8,
                          'Propane': 2.8, 'i-Butane': 4.8,
                          'n-Butane': 1.4,
                          'Nitrogen': 0, 'CO2': 0.61,
                          'H2S': (2 + 0.01 * p1 / 100)}
        pcp = 0
        for i in range(self.N_comp):
            if self.nm[i] in comp_plotn_std:
                pcp += comp_plotn_std[self.nm[i]] * conz1[i]

        p2 = 3900 # давление после дросселя, тдт, то, кПа
        T2, p2 = self.calc_drossel(p2)  # Температура газа после дросселя/тдт/то, С ; # Давление газа после дросселя/тдт/то, кПа

        C = 100  # Концентрация подаваемого метанола, %

        B = 0
        k = 0

        Pb = p2 / 100
        conz2 = [] # Состав газа после точки защиты
        for i in range(self.N_comp):
            conz2.append(self.gas[self.nm[i]])

        Z2 = self.H_y[2]  # К-т сжимаемости после

        Ggasstd = self.gas['Gm, kmole/h'] * 1000 * 24.465 / 1000  # ст. м3/ч при 25 С и 101.3 кПа

        if self.liq['n liq'] != 0:
            Mrliq2 = self.liq['Mr, kg/kmole']
            Gk2 = self.liq['G, kg/h'] / (Ggasstd / 1000)  # количество конденсата после точки защиты, кг/1000 м3
        else:
            Mrliq2 = 0
            Gk2 = 0

        'расчет температуры гидратообразования газа'

        if pcp == 0.55:
            B = 1.42
            k = 0.0193
        pass

        if 0.55 < pcp < 0.6:
            Bp = pcp - 0.55
            B = 1.42 - (Bp * 0.42 / 0.05)
            k = 0.0193 + (Bp * 0.0107 / 0.05)
        pass

        if pcp == 0.6:
            B = 1
            k = 0.03
        pass

        if 0.6 < pcp < 0.65:
            Bp = pcp - 0.6
            B = 1 - (Bp * 0.1 / 0.05)
            k = 0.03
        pass

        if pcp == 0.65:
            B = 0.9
            k = 0.03
        pass

        if 0.65 < pcp < 0.7:
            Bp = pcp - 0.65
            B = 0.9 - (Bp * 0.08 / 0.05)
            k = 0.03
        pass

        if pcp == 0.7:
            B = 0.82
            k = 0.03
        pass

        if 0.7 < pcp < 0.75:
            Bp = pcp - 0.7
            B = 0.82 - (Bp * 0.06 / 0.05)
            k = 0.03
        pass

        if pcp == 0.75:
            B = 0.76
            k = 0.03
        pass

        if 0.75 < pcp < 0.8:
            Bp = pcp - 0.75
            B = 0.76 - (Bp * 0.06 / 0.05)
            k = 0.03
        pass

        if pcp == 0.8:
            B = 0.7
            k = 0.03
        pass

        if 0.8 < pcp < 0.85:
            Bp = pcp - 0.8
            B = 0.7 - (Bp * 0.04 / 0.05)
            k = 0.03
        pass

        if pcp == 0.85:
            B = 0.66
            k = 0.03
        pass

        if 0.85 < pcp < 0.9:
            Bp = pcp - 0.85
            B = 0.66 - (Bp * 0.05 / 0.05)
            k = 0.03
        pass

        if pcp == 0.9:
            B = 0.61
            k = 0.03
        pass

        if 0.9 < pcp < 0.95:
            Bp = pcp - 0.9
            B = 0.61 - (Bp * 0.04 / 0.05)
            k = 0.03
        pass

        if pcp == 0.95:
            B = 0.57
            k = 0.03
        pass

        if 0.95 < pcp < 1:
            Bp = pcp - 0.95
            B = 0.57 - (Bp * 0.03 / 0.05)
            k = 0.03
        pass

        if pcp == 1:
            B = 0.54
            k = 0.03
        pass

        if pcp > 1:
            B = 0.54
            k = 0.03
        pass

        import math

        Thdr1 = (-(-0.035) - ((-0.035) ** 2 - 4 * (-0.035) * k * (math.log10(Pb) - round(B, 3))) ** (1 / 2)) / (
                2 * (-0.035) * k)
        Thdr2 = (-(-0.035) + ((-0.035) ** 2 - 4 * (-0.035) * k * (math.log10(Pb) - round(B, 3))) ** (1 / 2)) / (
                2 * (-0.035) * k)

        deltaT = 0

        if Thdr1 > Thdr2:
            Thdr = Thdr1
        else:
            Thdr = Thdr2

        if Thdr > 0 and T2 > 0:
            if T2 > Thdr:
                print("Температура гидратообразования, С = ", round(Thdr, 4))
                print("Гидраты не выпадают")
                deltaT = -1
            else:
                print("Температура гидратообразования, С = ", round(Thdr, 4))
                deltaT = Thdr - T2
        elif Thdr > 0 and T2 < 0:
            print("Температура гидратообразования, С = ", round(Thdr, 4))
            deltaT = Thdr - T2
        elif Thdr < 0 and T2 < 0:
            if T2 > Thdr:
                print("Температура гидратообразования, С = ", round(Thdr, 4))
                print("Гидраты не выпадают")
                deltaT = -1
            else:
                print("Температура гидратообразования, С = ", round(Thdr, 4))
                deltaT = Thdr - T2

        if deltaT > 0:
            e = 1
            X2 = 0
            while e > 0.01:
                if self.gas['Methane'] > 0.9:
                    A = 81 - 0.22 * X2 + 0.005 * X2 * (p2 / 100 - 7.5)
                else:
                    A = 81 - 0.33 * X2 + 0.01 * X2 * (p2 / 100 - 7.5)
                deT = -A * math.log((100 - X2) / (100 - 0.4378 * X2))
                e = abs(deltaT - deT)
                X2 += 0.001

            e1 = 1
            x = 0

            while e1 > 0.01:
                X22 = 3200 * x / (18 + 14 * x)
                x += 0.0001
                e1 = abs(X2 - X22)

            'расчет влагосодержания газа до и после точки подачи'

            lny1b1 = 2.4 - 530 / (T1 + 273.15)
            lny2b1 = 2.2 - 500 / (T1 + 273.15)

            y11 = math.exp(lny1b1 * (1 + (lny1b1 / lny2b1) * ((1 - x) / x)) ** (-2))
            y21 = math.exp(lny2b1 * (1 + (lny2b1 / lny1b1) * (x / (1 - x))) ** (-2))

            a11 = y11 * (1 - x)
            a21 = y21 * x

            lny1b2 = 2.4 - 530 / (T2 + 273.15)
            lny2b2 = 2.2 - 500 / (T2 + 273.15)

            y12 = math.exp(lny1b2 * (1 + (lny1b2 / lny2b2) * ((1 - x) / x)) ** (-2))
            y22 = math.exp(lny2b2 * (1 + (lny2b2 / lny1b2) * (x / (1 - x))) ** (-2))

            a12 = y12 * (1 - x)
            a22 = y22 * x

            comp_betta1 = {'Methane': math.exp(6.87 - 0.0093 * T1), 'Ethane': math.exp(4.649 - 0.0093 * T1),
                           'Propane': math.exp(7.665 - 0.00874 * T1), 'i-Butane': math.exp(7.91 - 0.00878 * T1),
                           'n-Butane': math.exp(7.98 - 0.0088 * T1),
                           'i-Pentane': math.exp(8.15 - 0.009 * T1), 'n-Pentane': math.exp(8.15 - 0.009 * T1),
                           'Nitrogen': math.exp(7.27 - 0.0012 * T1), 'CO2': math.exp(8.85 - 0.0117 * T1),
                           'H2S': math.exp(8.44 - 0.0091 * T1)}
            comp_betta2 = {'Methane': math.exp(6.87 - 0.0093 * T2), 'Ethane': math.exp(4.649 - 0.0093 * T2),
                           'Propane': math.exp(7.665 - 0.00874 * T2), 'i-Butane': math.exp(7.91 - 0.00878 * T2),
                           'n-Butane': math.exp(7.98 - 0.0088 * T2),
                           'i-Pentane': math.exp(8.15 - 0.009 * T2), 'n-Pentane': math.exp(8.15 - 0.009 * T2),
                           'Nitrogen': math.exp(7.27 - 0.0012 * T2), 'CO2': math.exp(8.85 - 0.0117 * T2),
                           'H2S': math.exp(8.44 - 0.0091 * T2)}
            comp_alpha1 = {'Methane': 0.725, 'Ethane': 0.6,
                           'Propane': 0.5, 'i-Butane': 0.4, 'n-Butane': 0.4,
                           'i-Pentane': 0.3, 'n-Pentane': 0.3, 'Nitrogen': 0.8, 'CO2': (0.568 - 0.0008 * T1),
                           'H2S': (0.56 - 0.0009 * T1)}
            comp_alpha2 = {'Methane': 0.725, 'Ethane': 0.6,
                           'Propane': 0.5, 'i-Butane': 0.4, 'n-Butane': 0.4,
                           'i-Pentane': 0.3, 'n-Pentane': 0.3, 'Nitrogen': 0.8, 'CO2': (0.568 - 0.0008 * T2),
                           'H2S': (0.56 - 0.0009 * T2)}

            alfa1 = []
            alfa2 = []
            betta1 = []
            betta2 = []
            for i in range(self.N_comp):
                if self.nm[i] in comp_betta1:
                    alfa1.append(comp_alpha1[self.nm[i]])
                    alfa2.append(comp_alpha2[self.nm[i]])
                    betta1.append(comp_betta1[self.nm[i]])
                    betta2.append(comp_betta2[self.nm[i]])
                else:
                    alfa1.append(0.0)
                    alfa2.append(0.0)
                    betta1.append(0.0)
                    betta2.append(0.0)

            a1cm = 0
            a2cm = 0
            b1cm = 0
            b2cm = 0

            for i in range(self.N_comp):
                a1cm += alfa1[i] * conz1[i]
                a2cm += alfa2[i] * conz2[i]
                b1cm += betta1[i] * conz1[i]
                b2cm += betta2[i] * conz2[i]

            p1 = p1 / 1000
            p2 = p2 / 1000
            W1 = 0.09984 * Z1 / p1 * math.exp(18.3036 - 3816.14 / (T1 + 273.15 - 46.13)) * math.exp(
                18 * p1 / (8.314 * (T1 + 273.15)) + (2 * b1cm * p1) / (Z1 * 8.314 * (T1 + 273.15) + a1cm * b1cm * p1))
            W2 = 0.09984 * Z2 / p2 * math.exp(18.3036 - 3816.14 / (T2 + 273.15 - 46.13)) * math.exp(
                18 * p2 / (8.314 * (T2 + 273.15)) + (2 * b2cm * p2) / (Z2 * 8.314 * (T2 + 273.15) + a2cm * b2cm * p2))

            'Расчет метанолосодержания'

            Ps1 = math.exp(18.5875 - 3626.55 / (T1 + 273.15 - 34.29)) * 1.3332 * 10 ** (-4)
            be1 = - math.exp(7.9154 - 0.01145 * (T1 + 273.15))

            Ps2 = math.exp(18.5875 - 3626.55 / (T2 + 273.15 - 34.29)) * 1.3332 * 10 ** (-4)
            be2 = - math.exp(7.9154 - 0.01145 * (T2 + 273.15))

            Q01 = 1331.31 * Ps1 / p1 * math.exp(-2 / 8.314 * p1 / (T1 + 273.15) * (be1 - 38.07 / 2))
            Q02 = 1331.31 * Ps2 / p2 * math.exp(-2 / 8.314 * p2 / (T2 + 273.15) * (be2 - 38.07 / 2))

            Q1 = a21 * Q01
            Q2 = a22 * Q02

            sigm1 = math.exp(-1.34 + 187.2 / C)
            sigm2 = math.exp(-1.34 + 187.2 / C)
            q1 = Gk1 * (2.6 - 0.016 * Mrliq1) * math.exp(12.2257 - 3903.6 / (T1 + 273.15) - sigm1)
            q2 = Gk2 * (2.6 - 0.016 * Mrliq2) * math.exp(12.2257 - 3903.6 / (T2 + 273.15) - sigm2)

            Gmethanol = (Gmeth1 * (X2 - X1) + (W1 - W2) * X2) / (C - X2) + (100 - X2) / (C - X2) * (
                        (Q2 - Q1) + (q2 - q1))

            if Gmethanol < 0:
                Gmethanol = 0
                print('Подача метанола не требуется')
            else:
                pass

            print("Расход метанола", round(Gmethanol, 4), "кг/1000 ст.м3 или", round(Gmethanol * Ggasstd / 1000, 4),
                  "кг/ч")
            return Gmethanol
        else:
            pass

    def comp_T_calc(self, p2, H2):
        def bisect(f, a, b, tol=1e-5):
            fa = f(a)
            fb = f(b)
            '''if fa * fb >= 0:
                raise ValueError("F")'''
            while (b - a) / 2 > tol:
                c = (a + b) / 2
                fc = f(c)
                if fc == 1:
                    return c
                if fa * fc < 0:
                    b = c
                    fb = fc
                else:
                    a = c
                    fa = fc
            return (a + b) / 2

        T = self.T
        def f(T):
            self.T = T
            self.p = p2
            self.calculation()
            H = self.H_calc()
            return abs(H2) - abs(H)

        root = bisect(f, -100, 100)
        print(root)
        return root

    def v_calc(self):
        # расчет молярного объема потока, м3/моль
        if self.e == 1:
            Zv = self.H_y[2]
            v = Zv * 8.314 * self.T_am / self.p # м3/моль
        else:
            Zl = self.H_x[2]
            v = Zl * 8.314 * self.T_am / self.p # м3/моль
        return v

    def ro_calc(self):
        # расчет плотности при рабочих условиях, г/см3
        if self.e == 1:
            Z = self.H_y[2]
        else:
            Z = self.H_x[2]
        ro = self.get_Mr() / (Z * self.R * self.T_am / self.p)
        return ro



    # производные уравнения Пенга-Робинсона

    def dpdV_T(self):
        # bar/(cm^3/mol)
        v = self.v_calc()
        if self.e == 1:
            B = self.Bv
            A = self.Av
        else:
            B = self.Bl
            A = self.Al
        dpdV_T = -self.R * self.T_am / (v - B)**2 + 2 * A * (v + B) / (v * (v + B) + B * (v - B))**2
        return dpdV_T
    
    def dadT_V(self, i):
        # cm^6 * bar/mol^2 * K
        dadT_V = - self.m[i] * self.Ap[i] / ((self.T_am * self.T_cr_am[i])**0.5 * (1+self.m[i]*(1-(self.T_am / self.T_cr_am[i])**0.5)))
        return dadT_V

    def dpdT_V(self, i):
        # bar/K
        v = self.v_calc() * 1000
        if self.e == 1:
            B = self.Bv
            A = self.Av
        else:
            B = self.Bl
            A = self.Al

        dpdT_V = self.R / (v - B) - self.dadT_V(i) / (v * (v + B) + B * (v - B))
        return dpdT_V
    
    def dadT_p(self, i):
        # bezrazm
        if self.e == 1:
            B = self.Bv
            A = self.Av
        else:
            B = self.Bl
            A = self.Al

        dadT_p = (((self.p/(self.R * self.T_am)**2) * (self.dadT_V(i) - 2*A/self.T_am)) * 10)/1000
        return dadT_p

    def dbdT_p(self):
        # bezrazm
        if self.e == 1:
            B = self.Bv
            A = self.Av
        else:
            B = self.Bl
            A = self.Al

        dbdT_p = - B * self.p_am / (self.R * self.T_am**2)
        return dbdT_p
    
    def dZdT_p(self, i):
        # bezrazm
        if self.e == 1:
            B = self.Bv
            A = self.Av
            Z = self.H_y[2]
        else:
            B = self.Bl
            A = self.Al
            Z = self.H_x[2]
        Av = 0.45723553 * self.alp[i] * (self.p_r[i] / self.T_r[i] ** 2)  # параметр А для газа в кубическом уравнении PR
        Bv = 0.077796074 * (self.p_r[i] / self.T_r[i])  # параметр B для газа в кубическом уравнении PR

        dZdT_p = (self.dadT_p(i) * (Bv - Z) + self.dbdT_p() * (6*Bv*Z + 2*Z - 3 * Bv**2 - 2*Bv + Av - Z**2)) / (3 * Z**2 + 2 * (Bv-1) * Z + (Av - 2*Bv - 3 * Bv**2))
        return dZdT_p
    
    def dvdT_p(self, i):
        # bezrazm
        if self.e == 1:
            B = self.Bv
            A = self.Av
            Z = self.H_y[2]
        else:
            B = self.Bl
            A = self.Al
            Z = self.H_x[2]

        dvdT_p = (self.R/self.p_am) * (self.T_am * self.dZdT_p(i) + Z)
        return dvdT_p
    
    def dTdp_V(self, i):
        # K/bar
        dTdp_V = 1/self.dpdT_V(i)
        return dTdp_V
    
    def d2adT2_v(self, i):
        # cm^6*bar/mol^2 * K^2
        d2adT2_v = self.m[i] * self.Ap[i] * (1+self.m[i]) * (self.T_cr_am[i] / self.T_am)**0.5 / (2*self.T_am * self.T_cr_am[i])
        return d2adT2_v


if __name__ == '__main__':

    comps = ['Methane', 'Ethane', 'Propane', 'H2O']
    concs = [0.5, 0.3, 0.15, 0.05]
    composition = {}
    for comp, conc in zip(comps, concs):
        composition.update({comp : conc})
    S1 = PengRobinson(composition=composition, pressure=101.3, temperature=20, mass_flow=100)

    #S1.phase_envelope()  #  функция для расчета фазовой кривой

    S1.calculation()  #  функция для расчета самого уравнения и составов фаз

    print("Энтальпия потока, кДж/кг", S1.H)
    print("Энтропия потока, кДж/кг*С", S1.S)
    print("Теплоемкость потока Cp, кДж/кг*С", S1.Cp)
    print("Теплоемкость потока Cv, кДж/кг*С", S1.Cv)
    print("Показатель адиабаты", S1.y)

    #S1.comp_T_calc(10000, -4495)  #  функция для расчета температуры на выходе из компрессора (первая цифры - новое давление, кПа; вторая цифра - новая энтальпия, кДж/кг)

    #S1.calc_drossel(3900)  #  функция для расчета дросселя и температуры на выходе из него (цифра - новое давление)

    #S1.Methanol_STO()  #  функция для расчета подачи метанола по СТО Газпром 3.1-3-010-2008

    if S1.gas['Gm, kmole/h'] != 0:
        print("Газовая фаза:", S1.gas)
    if S1.liq['Gm, kmole/h'] != 0:
        print("Жидкая фаза фаза:", S1.liq)
    if S1.water['Gm, kmole/h'] != 0:
        print("Водная фаза:", S1.water)

    # S1.phase_envelope()

# надо получить 7.768 кДж/кг С энтропия
# 3.205 кДж/кг С