// DewPointCpp.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include "Converters.hpp"
#include "DewPoint.hpp"

int main()
{
    std::vector<double> Tc{
        -82.45099487,  32.27800903, 96.74801025, 134.9460083,  152.04900513,
        187.24801025, 196.4500061, 234.74801025, 267.00802002, 295.44802246
    };
    std::vector<double> Pc{
        4640.68017578, 4883.85009766, 4256.66015625, 3647.62011719, 3796.62011719,
        3333.59008789, 3375.12011719, 3031.62011719, 2736.7800293,  2496.62011719
    };
    std::vector<double> Af{
        0.0114984,  0.0986,     0.1524,     0.18479,   0.20100001, 0.22224,
        0.25389001, 0.30070001, 0.34979001, 0.40180001
    };
    std::vector<double> Vc{
        0.0989999,  0.148,      0.2,        0.26300001, 0.25499001, 0.30799001,
        0.31099001, 0.368,      0.42598,    0.486
    };
    double Pressure = 101.325;
    std::vector<double> Yi(10, 0.1);    
    std::vector<double> Tr(Tc.size());
    double T = 273.15;
    
    UnitsConverter uc;
    DewPoint dp(Pressure, Yi, Tc, Pc, Af, Vc);
    double TDew = dp.Calculation();
    std::cout << uc.TemperatureUnits.RankineToCelcius(TDew);
    return 0;
};
