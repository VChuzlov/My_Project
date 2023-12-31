// DewPointCpp.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include "Converters.h"

int main()
{
    std::vector<double> Tc(10, 5.2);
    std::vector<double> Tr(Tc.size());
    double T = 273.15;
    
    ValuesConverter vc;
    Tr = vc.ReducedParam(T, Tc);

    for (int i = 0; i < Tc.size(); ++i)
    {
        std::cout << Tr[i] << " " << std::endl;
    }
    return 0;
};
