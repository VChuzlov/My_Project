#include <vector>
#include <math.h>
#include "DewPoint.hpp"

double DewPoint::CalculateAalpha(std::vector<double> mf,
    std::vector<std::vector<double>> kij, std::vector<double> ai,
    std::vector<double> alpha)
{
    double Result = 0.0;
    for (int i = 0; i < mf.size(); ++i)
    {
        for (int j = 0; j < mf.size(); ++j)
        {
            Result += (
                mf[i] * mf[j] * (1 - kij[i][j])
                * pow(ai[i] * alpha[i] * ai[j] * alpha[j], 0.5));
        }
    }
    return Result;
}

std::vector<double> DewPoint::CalculateM(std::vector<double> Af)
{
    std::vector<double> m(Af.size());
    return ;
}