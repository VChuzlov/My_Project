#include <vector>
#include <math.h>
#include <functional>
#include "DewPoint.hpp"
//#include "Functions.hpp" циклический импорт уже есть в DewPoint.hpp
#include "Converters.hpp"
#include <iostream>

double DewPoint::CalculateAalpha(
    const std::vector<double> &mf, 
    const std::vector<std::vector<double>> &kij, 
    const std::vector<double> &ai, const std::vector<double> &alpha)
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
    const std::vector<std::vector<double>> &kij, 
    const std::vector<double> &ap)
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

std::vector<double> DewPoint::CalculateAi(const std::vector<double> &tc,
    const std::vector<double> &pc)
{
    std::vector<double> Result(tc.size());
    for (unsigned int i = 0; i < tc.size(); ++i)
    {
        Result[i] = 0.45724 * pow(8.314 * tc[i], 2) / pc[i];
    }
    return Result;
}

double DewPoint::CalculateAl(const std::vector<double> &x,
    const std::vector<std::vector<double>> &ab)
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

std::vector<double> DewPoint::CalculateAlpha(const std::vector<double> &m,
    const std::vector<double> &tr)
{
    std::vector<double> Result(m.size());
    for (unsigned int i = 0; i < m.size(); ++i)
    {
        Result[i] = pow((1 + m[i] * (1 - pow(tr[i], 0.5))), 2);
    }
    return Result;
}

std::vector<double> DewPoint::CalculateAp(const std::vector<double> &alpha,
    const std::vector<double> &tr, const std::vector<double> &pr)
{
    std::vector<double> Result(tr.size());
    for (unsigned int i = 0; i < tr.size(); ++i)
    {
        Result[i] = 0.457235529 * alpha[i] * pr[i] / pow(tr[i], 2);
    }
    return Result;
}

double DewPoint::CalculateAv(const std::vector<double> &y,
    const std::vector<std::vector<double>> &ab)
{
    double Result = 0.0;
    for (unsigned int i = 0; i < y.size(); ++i)
    {
        for (unsigned int j = 0; j < y.size(); ++j)
        {
            Result += y[i] * y[j] * ab[i][j];
        }
    }
    return Result;
}

double DewPoint::CalculateBbl(
    const std::vector<double> &x, const std::vector<double> &bi)
{
    double Result = 0.0;
    for (unsigned int i = 0; i < x.size(); ++i)
    {
        Result += x[i] * bi[i];
    }
    return Result;
}

double DewPoint::CalculateBbv(
    const std::vector<double> &y, const std::vector<double> &bi)
{
    double Result = 0.0;
    for (unsigned int i = 0; i < y.size(); ++i)
    {
        Result += y[i] * bi[i];
    }
    return Result;
}

std::vector<double> DewPoint::CalculateBi(const std::vector<double> &tc,
    const std::vector<double> &pc)
{
    std::vector<double> Result(tc.size());
    for (unsigned int i = 0; i < tc.size(); ++i)
    {
        Result[i] = 0.07780 * 8.314 * tc[i] / pc[i];
    }
    return Result;
}

double DewPoint::CalculateBl(
    const std::vector<double> &x, const std::vector<double> &bp)
{
    double Result = 0.0;
    for (unsigned int i = 0; i < x.size(); ++i)
    {
        Result += x[i] * bp[i];
    }
    return Result;
}

std::vector<double> DewPoint::CalculateBp(const std::vector<double> &pr,
    const std::vector<double> &tr)
{
    std::vector<double> Result(pr.size());
    for (unsigned int i = 0; i < pr.size(); ++i)
    {
        Result[i] = 0.07796074 * pr[i] / tr[i];
    }
    return Result;
}

double DewPoint::CalculateBv(const std::vector<double> &y, 
    const std::vector<double> &bp)
{
    double Result = 0;
    for (unsigned int i = 0; i < y.size(); ++i)
    {
        Result += y[i] * bp[i];
    }
    return Result;
}

double DewPoint::CalculateD(
    const std::vector<double> &mf, const std::vector<double> &m,
    const std::vector<std::vector<double>> &kij, 
    const std::vector<double> &ai,
    const std::vector<double> &alpha, const std::vector<double> &tr)
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

std::vector<double> DewPoint::CalculateDi(const std::vector<double> &m,
    const std::vector<double> &ai, const std::vector<double> &alpha, 
    const std::vector<double> &tr)
{
    std::vector<double> Result(m.size());
    for (unsigned int i = 0; i < m.size(); ++i)
    {
        Result[i] = m[i] * ai[i] * alpha[i] * pow(tr[i] / alpha[i], 0.5);
    }
    return Result;
}

std::vector<double> DewPoint::CalculateFil(
    const std::vector<std::vector<double>> &ab, 
    const std::vector<double> &x, double zl,
    const std::vector<double> &bp, 
    double al, double bl)
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
    const std::vector<std::vector<double>> &ab, 
    const std::vector<double> &y, double zv,
    const std::vector<double> &bp, double av, double bv
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
    
    auto foo = [this](double t)
    {
        return this->ForInitialTValue(t);
    };

    Result = BrentsMethod(foo, 1e-5, 1200.0);
    return Result;
}

std::vector<std::vector<double>> DewPoint::CalculateKij(
    const std::vector<double> &vc, unsigned int n)
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

std::vector<double> DewPoint::CalculateXi(const std::vector<double> &ki)
{
    std::vector<double> Result(ki.size());
    double k = 1e-12;
    for (unsigned int i = 0; i < ki.size(); ++i)
    {
        ki[i] < 1e-12 ? k = k : k = ki[i];
        Result[i] = this->Yi[i] / k;
    }
    return Result;
}

double DewPoint::CalculateZl(
    double al, double bl, 
    std::function<std::vector<double> (double, double, double)> method)
{
    double Result = 0.0;
    std::vector<double> roots(3);
    roots = method(
        bl - 1.,
        al - 2. * bl - 3. * pow(bl, 2.),
        (-al + pow(bl, 2.) + bl) * bl
    );
    Result = this->SelectCubicEquationRoot(roots[0], roots[1], roots[2], Min);
    return Result;
}

double DewPoint::CalculateZv(
    double av, double bv,
    std::function<std::vector<double> (double, double, double)> method)
{
    double Result = 0.0;
    std::vector<double> roots(3);
    roots = method(
        bv - 1,
        av - 2 * bv - 3 * pow(bv, 2),
        (-av + pow(bv, 2) + bv) * bv
    );
    Result = this->SelectCubicEquationRoot(roots[0], roots[1], roots[2], Max);
    return Result;
}

double DewPoint::Calculation()
{
    double Result = 0.0;
    double T;
    std::vector<std::vector<double>> Kij;
    std::vector<double> m;
    UnitsConverter uc;
    unsigned int i = 0;
    double _t;

    T = this->CalculateInitialValueForT();
    Kij = this->CalculateKij(this->Vc);
    m = this->CalculateM(this->Af);
    this->PreCalculation(T, Kij, m);

    auto foo = [Kij, m, this](double t)
    {
        return this->InsideJob(t, Kij, m, this->XiNew);
    };

    while (!(this->Condition()))
    {
        i += 1;
        Result = BrentsMethod(
            foo,
            .95 * T,
            1.2 * T
        );
        _t = uc.TemperatureUnits.RankineToCelcius(Result);
        if (i > 10)
        {
            break;
        }
    }
    return Result;
}

bool DewPoint::Condition(double tol)
{
    unsigned int n = this->XiNew.size();
    double s = 0.0;

    for (unsigned int i = 0; i < n; ++i)
    {
        s += (pow(this->Xi[i] - this->XiNew[i], 2) / this->XiNew[i]) / n;
    }
    return s <= tol;
}

DewPoint::DewPoint(
    double pressure,
    const std::vector<double> &yi,
    const std::vector<double> &tc,
    const std::vector<double> &pc,
    const std::vector<double> &af,
    const std::vector<double> &volc)
{
    UnitsConverter uc;
    ValuesConverter vc;

    this->Pressure = uc.PressureUnits.kPaToPsi(pressure);
    this->Yi = yi;
    this->Af = af;
    this->Vc = volc;

    unsigned int size = yi.size();

    this->Tc.resize(size);
    this->Pc.resize(size);
    this->Tr.resize(size);
    this->Pr.resize(size);
    this->Xi.resize(size);
    this->XiNew.resize(size);

    for (unsigned int i = 0; i < size; ++i)
    {
        this->Tc[i] = uc.TemperatureUnits.CelciusToRankine(tc[i]);
        this->Pc[i] = uc.PressureUnits.kPaToPsi(pc[i]);
    }
    this->Pr = vc.ReducedParam(this->Pressure, this->Pc);
}

std::vector<double> DewPoint::EstimateKi(double t)
{
    std::vector<double> Result(this->Yi.size());
    for (unsigned int i = 0; i < this->Yi.size(); ++i)
    {
        Result[i] = exp(
            log(this->Pc[i] / this->Pressure)
            + log(10) * (7. / 3.) * (1. + this->Af[i])
            * (1. - this->Tc[i] / t)
        );
    }
    return Result;
}

double DewPoint::EstimateTFromXiAndTSati(
    const std::vector<double> &xi, 
    const std::vector<double> &tsati)
{
    double Result = 0.0;
    for (unsigned int i = 0; i < xi.size(); ++i)
    {
        Result += xi[i] * tsati[i];
    }
    return Result;
}

std::vector<double> DewPoint::EstimateTSati()
{
    std::vector<double> Result(this->Yi.size());
    for (unsigned int i = 0; i < Result.size(); ++i)
    {
        Result[i] = (
            this->Tc[i] / (1. - 3. * log(this->Pressure / this->Pc[i])
                           / (log(10) * (7. + 7. * this->Af[i]))
                          )
        );
    }
    return Result;
}

double DewPoint::ForInitialTValue(double t)
{
    std::vector<double> ki = this->EstimateKi(t);
    std::vector<double> xi = this->CalculateXi(ki);
    std::vector<double> tasti = this->EstimateTSati();
    double t_ = this->EstimateTFromXiAndTSati(xi, tasti);
    return t - t_;
}

double DewPoint::InsideJob(
    double t, const std::vector<std::vector<double>> &kij,
    const std::vector<double> &m, std::vector<double> &xi)
{
    double xSum = 0.0;
    for (unsigned int i = 0; i < xi.size(); ++i)
    {
        xSum += xi[i];
    }
    if (xSum != 1.0)
    {
        for (unsigned int i = 0; i < xi.size(); ++i)
        {
            xi[i] /= xSum;
        }
    }

    ValuesConverter vc;
    std::vector<double> Tr = vc.ReducedParam(t, this->Tc);
    std::vector<double> Alpha = this->CalculateAlpha(m, Tr);
    std::vector<double> Ap = this->CalculateAp(Alpha, Tr, this->Pr);
    std::vector<double> Bp = this->CalculateBp(this->Pr, Tr);
    std::vector<std::vector<double>> Ab = this->CalculateAb(kij, Ap);
    
    double Av = this->CalculateAv(this->Yi, Ab);
    double Bv = this->CalculateBv(this->Yi, Bp);
    double Zv = this->CalculateZv(Av, Bv);
    double Al = this->CalculateAl(xi, Ab);
    double Bl = this->CalculateBl(xi, Bp);
    double Zl = this->CalculateZl(Al, Bl);

    std::vector<double> Fiv = this->CalculateFiv(Ab, this->Yi, Zv, Bp, Av, Bv);
    std::vector<double> Fil = this->CalculateFil(Ab, xi, Zl, Bp, Al, Bl);

    double s = 0.0;
    std::vector<double> XiNew(m.size());
    for (unsigned int i = 0; i < m.size(); ++i)
    {
        XiNew[i] = this->Yi[i] * Fiv[i] / Fil[i];
        s += XiNew[i];
    }
    
    this->Xi = xi;
    this->XiNew = XiNew;
    
    return 1.0 - s;
}

void DewPoint::PreCalculation(
    double t, const std::vector<std::vector<double>> &kij,
    const std::vector<double> &m)
{
    ValuesConverter vc;
    std::vector<double> Tr = vc.ReducedParam(t, this->Tc);
    std::vector<double> Alpha = this->CalculateAlpha(m, Tr);
    std::vector<double> Ap = this->CalculateAp(Alpha, Tr, this->Pr);
    std::vector<double> Bp = this->CalculateBp(this->Pr, Tr);
    std::vector<std::vector<double>> Ab = this->CalculateAb(kij, Ap);
    std::vector<double> ki = this->EstimateKi(t);
    std::vector<double> xi = this->CalculateXi(ki);
    
    double Av = this->CalculateAv(this->Yi, Ab);
    double Bv = this->CalculateBv(this->Yi, Bp);
    double Zv = this->CalculateZv(Av, Bv);
    double Al = this->CalculateAl(xi, Ab);
    double Bl = this->CalculateBl(xi, Bp);
    double Zl = this->CalculateZl(Al, Bl);
    std::vector<double> Fiv = this->CalculateFiv(Ab, this->Yi, Zv, Bp, Av, Bv);
    std::vector<double> Fil = this->CalculateFil(Ab, xi, Zl, Bp, Al, Bl);

    std::vector<double> XiNew(m.size());
    for (unsigned int i = 0; i < m.size(); ++i)
    {
        XiNew[i] = this->Yi[i] * Fiv[i] / Fil[i];
    }

    this->Xi = xi;
    this->XiNew = XiNew;
}

double DewPoint::SelectCubicEquationRoot(
    double z1, double z2, double z3, 
    std::function<double (std::vector<double>)> foo)
{
    std::vector<double> roots{z1, z2, z3};
    return foo(roots);
}
