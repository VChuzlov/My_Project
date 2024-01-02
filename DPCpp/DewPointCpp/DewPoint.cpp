#include <vector>
#include <math.h>
#include "DewPoint.hpp"

double DewPoint::CalculateAalpha(std::vector<double> mf,
    std::vector<std::vector<double>> kij, std::vector<double> ai,
    std::vector<double> alpha)
{
    double Result = 0.0;
    for (unsigned int i = 0; i < mf.size(); ++i)
    {
        for (unsigned int j = 0; j < mf.size(); ++j)
        {
            Result += (
                mf[i] * mf[j] * (1 - kij[i][j])
                * pow(ai[i] * alpha[i] * ai[j] * alpha[j], 0.5));
        }
    }
    return Result;
}

std::vector<std::vector<double>> DewPoint::CalculateAb(
    std::vector<std::vector<double>> kij, std::vector<double> ap)
{
    unsigned int size = ap.size();
    std::vector<std::vector<double>> Result(size);
    for (unsigned int i = 0; i < size; ++i)
    {
        Result[i].resize(size);
        for (unsigned int j = 0; j < size; ++j)
        {
            Result[i][j] = (1 - kij[i][j]) * pow(ap[i] * ap[j], 0.5);
        }
    }
    return Result;
}

std::vector<double> DewPoint::CalculateAi(std::vector<double> tc,
    std::vector<double> pc)
{
    std::vector<double> Result(tc.size());
    for (unsigned int i = 0; i < tc.size(); ++i)
    {
        Result[i] = 0.45724 * pow(8.314 * tc[i], 2) / pc[i];
    }
    return Result;
}

double DewPoint::CalculateAl(std::vector<double> x,
    std::vector<std::vector<double>> ab)
{
    double Result = 0.0;
    for (unsigned int i = 0; i < x.size(); ++i)
    {
        for (unsigned int j = 0; j < x.size(); ++j)
        {
            Result += x[i] * x[j] * ab[i][j];
        }
    }
    return Result;
}

std::vector<double> DewPoint::CalculateM(std::vector<double> Af)
{
    std::vector<double> m(Af.size());
    return ;
}