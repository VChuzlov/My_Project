import math
class FlowPR:
    p = 7100  # давление, кПа
    T = 55  # температура, С
    N_comp = 6  # колиество компонентов
    R = 8.314   # Дж/(моль·K) газовая постоянная
    Name_comp = ['Метан', 'Этан', 'Пропан', 'Изобутан', 'Н-бутан', 'Азот']

    K = [0 for i in range(N_comp)]  # константы скорости
    p_cr = [4640.680, 4883.850, 4256.660, 3647.620, 3796.620, 3394.370]  # критические давления, кПа
    T_cr = [-82.45, 32.28, 96.75, 134.95, 152.05, -146.96]  # критические температуры, С
    w = [0.0115, 0.0986, 0.1524, 0.18479, 0.201, 0.04]  # ацентрические факторы
    conz = [0.7812, 0.03, 0.02, 0.0048, 0.001, 0.163]  # начальные конц,
    Mr = [16.0429000854492, 30.0699005126953, 44.0970001220703, 58.1240005493164, 58.1240005493164, 28.0130004882813]
    conz_x = [0 for i in range(N_comp)]  # концентрация в жид
    conz_y = [0 for i in range(N_comp)]  # концентрация в газ

    w_conz = [0 for i in range(N_comp)]  # массовая концентрация

    p_cr_am = [0 for i in range(N_comp)]  # критические давления, фунты на квадратный дюйм
    T_cr_am = [0 for i in range(N_comp)]  # критические температуры, градусы Рэнкина
    p_am = 0  # температура, градусы Рэнкина
    T_am = 0  # давление, фунты на квадратный дюйм
    p_r = [0 for i in range(N_comp)]  # приведенные давления
    T_r = [0 for i in range(N_comp)]  # приведенные температуры
    Ab = [[0 for i in range(6)] for j in range(N_comp)]
    kij = [[0, 0.00609, 0.01923, 0.02000, 0.02063, 0.02302],
           [0.00609, 0, -0.00297, 0.00031, 0.00019, 0.02464],
           [0.01923, -0.00297, 0, -0.00023, 0.00246, -0.01531],
           [0.02000, 0.00031, -0.00023, 0, -0.00109, 0.11068],
           [0.02063, 0.00019, 0.00246, -0.00109, 0, 0.00418],
           [0.02302, 0.02464, -0.01531, 0.11068, 0.00418, 0]
           ]  # бинарные коэффициенты

    a_Cp = [[2.36459,1.1429,0.395,0.1533,0.00854058,0.9827466],
            [-0.00426494,-0.0006472,0.00422818,0.00526958,0.00655398,0.00019428482],
            [0.0000169854,0.0000127293,0.000001189458,0.0000002181678,-0.00000332904,-0.0000000012473841],
            [-0.00000001489904,-0.00000001357264,-0.000000002668704,-0.000000002911584,0.000000000706584,-0.000000000014621936],
            [0.00000000000430448,0.00000000000441048,0.00000000000083968,0.00000000000118368,-3.19963E-14,2.0250665E-15]
           ]
    a_H = [[-12.98,-1.7675,39.4889,30.903,67.721,2.888634],
           [2.36459,1.1429,0.395,0.1533,0.00854058,0.9827466],
           [-0.00213247,-0.0003236,0.00211409,0.00263479,0.00327699,0.00009714241],
           [0.0000056618,0.0000042431,0.000000396486,0.0000000727226,-0.00000110968,-0.0000000004157947],
           [-0.00000000372476,-0.00000000339316,-0.000000000667176,-0.000000000727896,0.000000000176646,-0.000000000003655484],
           [0.000000000000860896,0.000000000000882096,0.000000000000167936,0.000000000000236736,-6.39926E-15,4.050133E-16]
          ]
    Bp = [0 for i in range(N_comp)]
    Ap = [0 for i in range(N_comp)]
    Av = 0
    Bv = 0
    Al = 0
    Bl = 0
    e = 0  # доля отгона
    Zv = 0
    Zl = 0
    Fiv = [0 for i in range(N_comp)]
    Fil = [0 for i in range(N_comp)]
    OldK = [0 for i in range(N_comp)]

    def recalc_P_T(self):
        '''Пересчет температуры в град. Рэнкина, давления в фунты на квадратный дюйм'''
        self.p_am = self.p * 0.145
        self.T_am = self.T * 1.8 + 491.67
        for i in range(self.N_comp):
            self.p_cr_am[i] = self.p_cr[i] * 0.145
            self.T_cr_am[i] = self.T_cr[i] * 1.8 + 491.67
            self.p_r[i] = self.p_am / self.p_cr_am[i]
            self.T_r[i] = self.T_am / self.T_cr_am[i]

    def calcK_Wilson(self):
        '''Расчет констант равновесия по уравнению Вильсона для начального приближения'''
        for i in range(self.N_comp):
            self.K[i] = self.p_cr_am[i] / self.p_am * math.exp(
                5.37 * (1 + self.w[i]) * (1 - self.T_cr_am[i] / self.T_am))

    def calc_Ab_Bp(self):
        '''Расчет вириальных коэффициентов'''
        Alfa = [0 for i in range(self.N_comp)]
        m = [0 for i in range(self.N_comp)]
        for i in range(self.N_comp):
            self.Bp[i] = 0.077796074 * self.p_r[i] / self.T_r[i]
            m[i] = 0.37464 + 1.54226 * self.w[i] - 0.26992 * self.w[i] * self.w[i]
            Alfa[i] = math.pow(1 + m[i] * (1 - math.sqrt(self.T_r[i])), 2)
            self.Ap[i] = 0.457235529 * Alfa[i] * self.p_r[i] / math.pow(self.T_r[i], 2)
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

    def FazeCalc(self):
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

    def Rashford_Rice(self, _fc):  # расчет мольной доли отгона
        def fn(_e):
            s = 0;
            for i in range(self.N_comp):
                if self.K[i] > 0:
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
                _e = (a + b) / 2;
                if fn(a) * fn(_e) > 0:
                    a = _e
                else:
                    b = _e
                _e = (a + b) / 2

        for i in range(self.N_comp):
            self.conz_x[i] = self.conz[i] / (1 + _e * (self.K[i] - 1));
            self.conz_y[i] = self.K[i] * self.conz_x[i]

        self.e = _e;

    def calc_x_y(self):

        fc = self.FazeCalc()
        e_old = self.e
        if fc in [0, 3, 4]:
            self.Rashford_Rice(fc)
            if self.e >= 1:
                self.e = 1
                for i in range(self.N_comp):
                    self.conz_x[i] = 0
                    self.conz_y[i] = self.conz[i]

            if self.e <= 0:
                self.e = 0;
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

    def calc_Fiv_Fil(self):
        for i in range(self.N_comp):
            sum_x = 0
            sum_y = 0
            for j in range(self.N_comp):
                sum_x += self.Ab[i][j] * self.conz_x[j]
                sum_y += self.Ab[i][j] * self.conz_y[j]

            self.Fiv[i] = math.exp((self.Zv - 1) * self.Bp[i] / self.Bv - math.log(self.Zv - self.Bv) - \
                                   self.Av / (math.sqrt(2) * self.Bv) * (sum_y / self.Av - \
                                                                         self.Bp[i] / self.Bv / 2) * math.log(
                (self.Zv + (1 + math.sqrt(2)) * self.Bv) / (self.Zv - (-1 + math.sqrt(2)) * self.Bv)))

            self.Fil[i] = math.exp((self.Zl - 1) * self.Bp[i] / self.Bl - math.log(self.Zl - self.Bl) - \
                                   self.Al / (math.sqrt(2) * self.Bl) * (sum_x / self.Al - \
                                                                         self.Bp[
                                                                             i] / self.Bl / 2) * math.log(
                (self.Zl + (1 + math.sqrt(2)) * self.Bl) / (self.Zl - (-1 + math.sqrt(2)) * self.Bl)))

    def memoryK(self):
        for i in range(self.N_comp):
            self.OldK[i] = self.K[i]

    def calcK_PR(self):
        for i in range(self.N_comp):
            self.K[i] = self.Fil[i] / self.Fiv[i]

    def calc_pogr(self):
        s = 0
        for i in range(self.N_comp):
            s += abs(self.OldK[i] - self.K[i])
        return s

    def pred_calculation(self):
        self.recalc_P_T()
        self.calcK_Wilson()
        self.memoryK()

    def _calculation(self):
        def func(Z1, Z2, Z3):
            if Z1 != 'null':
                return Z1
            if Z2 != 'null':
                return Z2
            if Z3 != 'null':
                return Z3

        self.calc_x_y()
        # теперь мы знаем конц x и y, далее опять по статье
        self.calc_Ab_Bp()
        self.calc_Av_Bv()
        # куб уравнение
        Z1, Z2, Z3 = self.calcCubeExp(self.Bv - 1, self.Av - 2 * self.Bv - 3 * math.pow(self.Bv, 2),
                                      (-self.Av + math.pow(self.Bv, 2) + self.Bv) * self.Bv)
        self.Zv = func(Z1, Z2, Z3)
        if (self.Zv==1):
            return False
        self.calc_Al_Bl()
        # куб уравнение
        Z1, Z2, Z3 = self.calcCubeExp(self.Bl - 1, self.Al - 2 * self.Bl - 3 * math.pow(self.Bl, 2),
                                      (-self.Al + math.pow(self.Bl, 2) + self.Bl) * self.Bl)
        self.Zl = func(Z1, Z2, Z3)
        if (self.Zl==1):
            return False
        self.calc_Fiv_Fil()
        self.calcK_PR()
        return True

    def calculation(self):
        self.pred_calculation()
        self._calculation()
        while self.calc_pogr() > 0.0001:
            self.memoryK()
            if not self._calculation():
                break
        print(self.K)
        print(self.conz_x)
        print(self.conz_y)

    def calc_k_J_T(self):
        '''Расчет коэффициента Джоуля-Томсона'''
        self.calculation()
        z1 = self.Zv
        t0  = self.T
        alfa = 0.0001
        self.T = self.T+alfa
        self.calculation()
        z2 = self.Zv
        #http://vesti-gas.ru/sites/default/files/attachments/vgn-1-42-2020-023-031.pdf
        # https://neftegas.info/upload/iblock/4d9/4d9537ab49679682f7211efb2ae23517.pdf
        mu = (t0+273.15)**2*self.R/(self.p/1000*self.get_Cp()*1000)*(z2-z1)/alfa +(t0+273.15)*self.R*z1/(self.p/1000*self.get_Cp()*1000)
        print(mu)

    def get_Cp(self):
        '''Функция для расчета '''
        self.calc_conz_x_to_w()
        Cv = 0
        for i in range(self.N_comp):
            Cv += self.w_conz[i]*(self.a_Cp[0][i]+self.a_Cp[1][i]*(self.T+273.15)+self.a_Cp[2][i]*(self.T+273.15)**2+
                                self.a_Cp[3][i] * (self.T + 273.15)**3 + self.a_Cp[4][i]*(self.T+273.15)**4)
        Mr_s = self.get_Mr()
        Cp = (Cv*Mr_s+self.R)/Mr_s
        return Cp

    def get_Mr(self):
        '''Функция для расчета средней молекулярной массы'''
        Mr_s = 0
        for i in range(self.N_comp):
            Mr_s += self.conz[i] * self.Mr[i]
        return Mr_s

    def get_H(self):
        '''Функция для расчета энтальпии'''
        self.calc_conz_x_to_w()
        H = 0
        for i in range(self.N_comp):
            H += self.w_conz[i] * (self.a_H[0][i] + self.a_H[1][i] * (self.T + 273.15) + self.a_H[2][i] * (
                        self.T + 273.15) ** 2 +
                                  self.a_H[3][i] * (self.T + 273.15) ** 3 + self.a_H[4][i] * (self.T + 273.15) ** 4
                                 + self.a_H[5][i] * (self.T + 273.15) ** 5)
        return H

    def calc_conz_x_to_w(self):
        '''Пересчет из мольнхых долей в массовые'''
        Mr_s = self.get_Mr()
        for i in range(self.N_comp):
            self.w_conz[i] = self.Mr[i]/Mr_s*self.conz[i]

    def calc_mathanol(self):
        # Tp = float(input("Введите температуру рабочего потока, C: "))'
        # u = float(input("Введите коэффициент эжекции (расход низконапорного газа к рабочему газу): "))'
        # Tn = float(input("Введите температуру низконапорного потока, C: "))'
        # Pp = float(input("Введите давление рабочего газа, МПа: "))'
        # Pc = float(input("Введите давление газа на выходе из эжектора, МПа: "))'
        # Cp = float(input("Введите среднюю изобарную теплоемкость газа, кДж/кг*К: "))'
        Tp = 35
        u = 1
        Tn = 22
        Pp = 8
        Pc = 6
        'Cp = 2.23'
        E1 = 0.98 * 10 ** 6
        E2 = 1.5
        Tp_ = Tp + 273
        Tn_ = Tn + 273

        Cp = 1.696 + 1.838 * 10 ** (-3) * Tp_ + 1.96 * (10 ** 6) * (Pp - 1) / Tp_ ** 3
        D = (1 / (Cp * 1000)) * (E1 / Tp ** 2 - E2)
        print("коэффициент Джоуля-Томпсона равен", round(D, 3), "К/МПа")

        Tc = (Tp_ + u * Tn_ - D * (Pp - Pc)) / (1 + u)
        Tc1 = Tc - 273
        print("температура на выходе из эжектора равна", round(Tc1))

        Db = 0.5
        Dm = 0.2
        V = 150000  # м3/ч
        g = 9.81
        Pl = 0.72  # кг/м3
        Pn = 2
        u1 = (V/3600)/(0.785*Db**2)
        u2 = (V/3600)/(0.785*Dm**2)
        DeltaP = ((Pl*(u2**2-u1**2))/2)*10**(-3)
        print(DeltaP, "кПа")

        Pej = Pp-(DeltaP*10**(-3))

        if Pej<Pn:
           print("эжектор работает")
        else:
           print("эжектор не работает")


if __name__ == '__main__':
    F = FlowPR()
    '''F.calculation()
    ee = []
    ee1 = []
    for i in range(-70, -52):
        F.T = i
        F.calculation()
        ee.append(i)
        ee1.append(F.e)
    print(ee)
    print(ee1)'''
    F.calc_k_J_T()
    #F.calc_mathanol()


