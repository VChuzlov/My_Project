#include <vector>
#include <math.h>
#include <functional>
#include "DewPoint.hpp"
#include "Functions.hpp"

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

std::vector<double> DewPoint::CalculateAlpha(std::vector<double> m,
    std::vector<double> tr)
{
    std::vector<double> Result(m.size());
    for (unsigned int i = 0; i < m.size(); ++i)
    {
        Result[i] = pow((1 + m[i] * (1 - pow(tr[i], 0.5))), 2);
    }
    return Result;
}

std::vector<double> DewPoint::CalculateAp(std::vector<double> alpha,
    std::vector<double> tr, std::vector<double> pr)
{
    std::vector<double> Result(tr.size());
    for (unsigned int i = 0; i < tr.size(); ++i)
    {
        Result[i] = 0.457235529 * alpha[i] * pr[i] / pow(tr[i], 2);
    }
    return Result;
}

double DewPoint::CalculateAv(std::vector<double> y,
    std::vector<std::vector<double>> ab)
{
    double Result;
    for (unsigned int i = 0; i < y.size(); ++i)
    {
        for (unsigned int j = 0; j < y.size(); ++j)
        {
            Result += y[i] * y[j] * ab[i][j];
        }
    }
    return Result;
}

double DewPoint::CalculateBbl(std::vector<double> x, std::vector<double> bi)
{
    double Result = 0.0;
    for (unsigned int i = 0; i < x.size(); ++i)
    {
        Result += x[i] * bi[i];
    }
    return Result;
}

double DewPoint::CalculateBbv(std::vector<double> y, std::vector<double> bi)
{
    double Result = 0.0;
    for (unsigned int i = 0; i < y.size(); ++i)
    {
        Result += y[i] * bi[i];
    }
    return Result;
}

std::vector<double> DewPoint::CalculateBi(std::vector<double> tc,
    std::vector<double> pc)
{
    std::vector<double> Result(tc.size());
    for (unsigned int i = 0; i < tc.size(); ++i)
    {
        Result[i] = 0.07780 * 8.314 * tc[i] / pc[i];
    }
    return Result;
}

double DewPoint::CalculateBl(std::vector<double> x, std::vector<double> bp)
{
    double Result = 0.0;
    for (unsigned int i = 0; i < x.size(); ++i)
    {
        Result += x[i] * bp[i];
    }
    return Result;
}

std::vector<double> DewPoint::CalculateBp(std::vector<double> pr,
    std::vector<double> tr)
{
    std::vector<double> Result(pr.size());
    for (unsigned int i = 0; i < pr.size(); ++i)
    {
        Result[i] = 0.07796074 * pr[i] / tr[i];
    }
    return Result;
}

double DewPoint::CalculateBv(std::vector<double> y, std::vector<double> bp)
{
    double Result = 0;
    for (unsigned int i = 0; i < y.size(); ++i)
    {
        Result += y[i] * bp[i];
    }
    return Result;
}

double DewPoint::CalculateD(std::vector<double> mf, std::vector<double> m,
    std::vector<std::vector<double>> kij, std::vector<double> ai,
    std::vector<double> alpha, std::vector<double> tr)
{
    double Result = 0.0;
    for (unsigned int i = 0; i < m.size(); ++i)
    {
        for (unsigned int j = 0; j < m.size(); ++j)
        {
            Result += (
                mf[i] * mf[j] * m[j] * (1 - kij[i][j])
                * pow(ai[i] * alpha[i] * alpha[j] * tr[j], 0.5) 
            );
        }
    }
    return Result;
}

std::vector<double> DewPoint::CalculateDi(std::vector<double> m,
    std::vector<double> ai, std::vector<double> alpha, std::vector<double> tr)
{
    std::vector<double> Result(m.size());
    for (unsigned int i = 0; i < m.size(); ++i)
    {
        Result[i] = m[i] * ai[i] * alpha[i] * pow(tr[i] / alpha[i], 0.5);
    }
    return Result;
}

std::vector<double> DewPoint::CalculateFil(
    std::vector<std::vector<double>> ab, std::vector<double> x, double zl,
    std::vector<double> bp, double al, double bl)
{
    std::vector<double> Result(x.size());
    double s;
    for (unsigned int i = 0; i < x.size(); ++i)
    {
        s = 0.0;
        for (unsigned int j = 0; j < x.size(); ++j)
        {
            s += ab[i][j] * x[j];
        }
        Result[i] = exp(
            (zl - 1) * bp[i] / bl - log(zl - bl)
            - al / (2 * pow(2, 0.5) * bl)
            * (2 * s / al - bp[i] / bl)
            * log((zl + (1 + pow(2, 0.5)) * bl)
                    / (zl - (-1 + pow(2, 0.5)) * bl))
        );
    }
    return Result;
}

std::vector<double> DewPoint::CalculateFiv(
    std::vector<std::vector<double>> ab, std::vector<double> y, double zv,
    std::vector<double> bp, double av, double bv
)
{
    std::vector<double> Result(y.size());
    double s;
    for (unsigned int i = 0; i < y.size(); ++i)
    {
        s = 0.0;
        for (unsigned int j = 0; j < y.size(); ++j)
        {
            s += ab[i][j] * y[j];
        }
        Result[i] = exp(
            (zv - 1) * bp[i] / bv - log(zv - bv)
            - av / (2 * pow(2, 0.5) * bv)
            * (2 * s / av - bp[i] / bv)
            * log((zv + (1 + pow(2, 0.5)) * bv)
                    / (zv - (-1 + pow(2, 0.5)) * bv))
        );
    }
    return Result;
}

double DewPoint::CalculateInitialValueForT()
{
    double Result = 0.0; 
    Result = BrentsMethod(ForinitialTValue, 1e-5, 1000.0);
    return Result;
}

std::vector<std::vector<double>> DewPoint::CalculateKij(
    std::vector<double> vc, unsigned int n)
{
    unsigned int size = vc.size();
    std::vector<double> VcR3(size);
    std::vector<std::vector<double>> Numerator(size);
    std::vector<std::vector<double>> Denominator(size);
    std::vector<std::vector<double>> Result(size);
    for (unsigned int i = 0; i < size; ++i)
    {
        Result[i].resize(size);
        Numerator[i].resize(size);
        Denominator[i].resize(size);
        VcR3[i] = pow(vc[i], 1. / 3.);
    }
    for (unsigned int i = 0; i < size; ++i)
    {
        for (unsigned int j = 0; j < size; ++j)
        {
            Numerator[i][j] = pow(VcR3[i] * VcR3[j], 0.5);
            Denominator[i][j] = (VcR3[i] + VcR3[j]) / 2.;
            Result[i][j] = 1 - pow(Numerator[i][j] / Denominator[i][j], n);
        }
    }
    return Result;
}

std::vector<double> DewPoint::CalculateM(std::vector<double> af)
{
    std::vector<double> Result(af.size());
    for (unsigned int i = 0; i < af.size(); ++i)
    {
        if (af[i] <= 0.49)
        {
            Result[i] = 0.3796 + 1.54226 * af[i] - 0.26992 * pow(af[i], 2.);
        }
        else
        {
            Result[i] = (0.379642 + 1.48503 * af[i] - 0.1644 * pow(af[i], 2.)
                         + 0.016667 * pow(af[i], 3.));
        }
    }
    return Result;
}