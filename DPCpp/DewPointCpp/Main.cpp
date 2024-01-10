// DewPointCpp.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <chrono>
#include "Converters.hpp"
#include "DewPoint.hpp"

int main()
{
    std::vector<double> Tc{
        -82.451,   32.278,   96.748,  134.946,  152.049,
        187.248,  196.450,  234.748,  267.008,  295.448,
        321.448,  344.448,  365.149,  385.149,  402.649,
        420.850,  433.850,  443.850,  460.220,  472.110,
        482.779,  494.850,  504.850,  513.850,  522.850,
        530.850,  538.850,  545.850,  552.850,  558.850
    };
    std::vector<double> Pc{
        4640.680, 4883.850, 4256.660, 3647.620, 3796.620,
        3333.590, 3375.120, 3031.620, 2736.780, 2496.620,
        2300.070, 2107.550, 1964.930, 1829.920, 1723.530,
        1620.180, 1516.810, 1420.560, 1316.900, 1213.470,
        1116.950, 1160.000, 1110.000, 1060.000, 1020.000,
        980.000,  950.000,  910.000,  883.000,  850.000
    };
    std::vector<double> Af{
        0.011,    0.099,    0.152,    0.185,    0.201,
        0.222,    0.254,    0.301,    0.350,    0.402,
        0.445,    0.488,    0.535,    0.562,    0.623,
        0.679,    0.706,    0.765,    0.770,    0.800,
        0.827,    0.907,    0.942,    0.972,    1.026,
        1.071,    1.105,    1.154,    1.214,    1.238
    };
    std::vector<double> Vc{
        0.099,    0.148,    0.200,    0.263,    0.255,
        0.308,    0.311,    0.368,    0.426,    0.486,
        0.543,    0.602,    0.660,    0.713,    0.780,
        0.830,    0.880,    0.940,    1.006,    1.070,
        1.130,    1.190,    1.245,    1.301,    1.356,
        1.411,    1.465,    1.518,    1.571,    1.622
    };
    double Pressure = 101.325;
    std::vector<double> Yi(30, 0.0333);    
    std::vector<double> Tr(Tc.size());
    double T = 273.15;
    
    UnitsConverter uc;
    DewPoint dp(Pressure, Yi, Tc, Pc, Af, Vc);

    auto start = std::chrono::high_resolution_clock::now();
    double TDew = dp.Calculation();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    std::cout << uc.TemperatureUnits.RankineToCelcius(TDew) << std::endl;
    std::cout << "Elapsed time: " << duration.count() / 1e6 << " seconds" << std::endl;
    return 0;
};
