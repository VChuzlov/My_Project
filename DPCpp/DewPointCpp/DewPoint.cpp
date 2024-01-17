#include <math.h>
#include <functional>
#include "DewPoint.hpp"
//#include "Functions.hpp" циклический импорт уже есть в DewPoint.hpp
#include "Converters.hpp"
#include <iostream>

void DewPoint::CalculateAalpha(const std::vector<double> &mf)
{
    this->Aalpha = 0.0;
    for (size_t i = 0; i < mf.size(); ++i)
    {
        for (size_t j = 0; j < mf.size(); ++j)
        {
            this->Aalpha += (
                mf[i] * mf[j] * (1 - this->Kij[i][j])
                * pow(this->Ai[i] * this->Alpha[i] 
                        * this->Ai[j] * this->Alpha[j], 0.5));
        }
    }
}

void DewPoint::CalculateAb(
    const std::vector<std::vector<double>> &kij, 
    const std::vector<double> &ap)
{
    size_t size = ap.size();
    for (size_t i = 0; i < size; ++i)
    {
        for (size_t j = 0; j < size; ++j)
        {
            this->Ab[i][j] = (1 - kij[i][j]) * pow(ap[i] * ap[j], 0.5);
        }
    }
}

void DewPoint::CalculateAi(const std::vector<double> &tc,
    const std::vector<double> &pc)
{
    for (size_t i = 0; i < tc.size(); ++i)
    {
        this->Ai[i] = 0.45724 * pow(8.314 * tc[i], 2) / pc[i];
    }
}

void DewPoint::CalculateAl(const std::vector<double> &x)
{
    this->Al = 0.0;
    for (size_t i = 0; i < x.size(); ++i)
    {
        for (size_t j = 0; j < x.size(); ++j)
        {
            this->Al += x[i] * x[j] * this->Ab[i][j];
        }
    }
}

void DewPoint::CalculateAlpha(const std::vector<double> &m,
    const std::vector<double> &tr)
{
    for (size_t i = 0; i < m.size(); ++i)
    {
        this->Alpha[i] = pow((1 + m[i] * (1 - pow(tr[i], 0.5))), 2);
    }
}

void DewPoint::CalculateAp(const std::vector<double> &alpha,
    const std::vector<double> &tr, const std::vector<double> &pr)
{
    for (size_t i = 0; i < tr.size(); ++i)
    {
        this->Ap[i] = 0.457235529 * alpha[i] * pr[i] / pow(tr[i], 2);
    }
}

void DewPoint::CalculateAv()
{
    this->Av = 0.0;
    for (size_t i = 0; i < this->Yi.size(); ++i)
    {
        for (size_t j = 0; j < this->Yi.size(); ++j)
        {
            this->Av += this->Yi[i] * this->Yi[j] * this->Ab[i][j];
        }
    }
}

void DewPoint::CalculateBbl(
    const std::vector<double> &x)
{
    this->Bbl = 0.0;
    for (size_t i = 0; i < x.size(); ++i)
    {
        this->Bbl += x[i] * this->Bi[i];
    }
}

void DewPoint::CalculateBbv()
{
    this->Bbv = 0.0;
    for (size_t i = 0; i < this->Yi.size(); ++i)
    {
        this->Bbv += this->Yi[i] * this->Bi[i];
    }
}

void DewPoint::CalculateBi(const std::vector<double> &tc,
    const std::vector<double> &pc)
{
    for (size_t i = 0; i < tc.size(); ++i)
    {
        this->Bi[i] = 0.07780 * 8.314 * tc[i] / pc[i];
    }
}

void DewPoint::CalculateBl(const std::vector<double> &x)
{
    this->Bl = 0.0;
    for (size_t i = 0; i < x.size(); ++i)
    {
        this->Bl += x[i] * this->Bp[i];
    }
}

void DewPoint::CalculateBp(const std::vector<double> &pr,
    const std::vector<double> &tr)
{
    for (size_t i = 0; i < pr.size(); ++i)
    {
        this->Bp[i] = 0.07796074 * pr[i] / tr[i];
    }
}

void DewPoint::CalculateBv()
{
    this->Bv = 0;
    for (size_t i = 0; i < this->Yi.size(); ++i)
    {
        this->Bv += this->Yi[i] * this->Bp[i];
    }
}

void DewPoint::CalculateD(
    const std::vector<double> &mf, const std::vector<double> &tr)
{
    this->D = 0.0;
    for (size_t i = 0; i < mf.size(); ++i)
    {
        for (size_t j = 0; j < mf.size(); ++j)
        {
            this->D += (
                mf[i] * mf[j] * this->M[j] * (1 - this->Kij[i][j])
                * pow(this->Ai[i] * this->Alpha[i] 
                        * this->Alpha[j] * tr[j], 0.5) 
            );
        }
    }
}

void DewPoint::CalculateDi(const std::vector<double> &tr)
{
    for (size_t i = 0; i < tr.size(); ++i)
    {
        this->Di[i] = (
            this->M[i] * this->Ai[i] * this->Alpha[i] 
            * pow(tr[i] / this->Alpha[i], 0.5)
        );
    }
}

std::vector<double> DewPoint::CalculateFil(
    const std::vector<std::vector<double>> &ab, 
    const std::vector<double> &x, const double &zl,
    const std::vector<double> &bp, 
    const double &al, const double &bl)
{
    std::vector<double> Result(x.size());
    double s;
    for (size_t i = 0; i < x.size(); ++i)
    {
        s = 0.0;
        for (size_t j = 0; j < x.size(); ++j)
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
    const std::vector<double> &y, const double &zv,
    const std::vector<double> &bp, 
    const double &av, const double &bv
)
{
    std::vector<double> Result(y.size());
    double s;
    for (size_t i = 0; i < y.size(); ++i)
    {
        s = 0.0;
        for (size_t j = 0; j < y.size(); ++j)
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

void DewPoint::CalculateKij(const std::vector<double> &vc, int n)
{
    size_t size = vc.size();
    std::vector<double> VcR3(size);
    
    for (size_t i = 0; i < size; ++i)
    {
        VcR3[i] = pow(vc[i], 1. / 3.);
    }
    for (size_t i = 0; i < size; ++i)
    {
        for (size_t j = 0; j < size; ++j)
        {
            if (n == 1)
            {
                this->Kij[i][j] = (
                    1 - pow(VcR3[i] * VcR3[j], 0.5) 
                    / ((VcR3[i] + VcR3[j]) / 2.));
            }
            else
            {
                this->Kij[i][j] = 1 - pow(
                    pow(VcR3[i] * VcR3[j], 0.5) 
                    / ((VcR3[i] + VcR3[j]) / 2.), n);
            }
            
        }
    };
}

void DewPoint::CalculateM(
    const std::vector<double> &af)
{
    for (size_t i = 0; i < af.size(); ++i)
    {
        if (af[i] <= 0.49)
        {
            this->M[i] = 0.3796 + 1.54226 * af[i] - 0.26992 * pow(af[i], 2.);
        }
        else
        {
            this->M[i] = (0.379642 + 1.48503 * af[i] - 0.1644 * pow(af[i], 2.)
                         + 0.016667 * pow(af[i], 3.));
        }
    }
}

std::vector<double> DewPoint::CalculateXi(const std::vector<double> &ki)
{
    std::vector<double> Result(ki.size());
    double k = 1e-12;
    for (size_t i = 0; i < ki.size(); ++i)
    {
        ki[i] < 1e-12 ? k = k : k = ki[i];
        Result[i] = this->Yi[i] / k;
    }
    return Result;
}

double DewPoint::CalculateZl(
    const double &al, const double &bl, 
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

void DewPoint::CalculateZv(
    std::function<std::vector<double> (double, double, double)> method)
{
    this->Zv = 0.0;
    std::vector<double> roots(3);
    roots = method(
        this->Bv - 1,
        this->Av - 2 * this->Bv - 3 * pow(this->Bv, 2),
        (-this->Av + pow(this->Bv, 2) + this->Bv) * this->Bv
    );
    this->Zv = this->SelectCubicEquationRoot(roots[0], roots[1], roots[2], Max);
}

double DewPoint::Calculation()
{
    double Result = 0.0;
    double T;
    const size_t size = this->Vc.size();
    
    UnitsConverter uc;
    size_t i = 0;
    double _t;

    T = this->CalculateInitialValueForT();
    this->CalculateKij(this->Vc);
    this->CalculateM(this->Af);
    this->PreCalculation(T, this->Kij, this->M);

    auto foo = [this](double t)
    {
        return this->InsideJob(t, this->Kij, this->M, this->XiNew);
    };

    while (!(this->Condition()))
    {
        i += 1;
        Result = BrentsMethod(
            foo,
            .95 * T,
            1.2 * T
        );
        if (i > 1000)
        {
            break;
        }
    }
    return Result;
}

bool DewPoint::Condition(double tol)
{
    size_t n = this->XiNew.size();
    double s = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        s += (pow(this->Xi[i] - this->XiNew[i], 2) / this->XiNew[i]) / n;
    }
    return s <= tol;
}

DewPoint::DewPoint(
    const double &pressure,
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

    size_t size = yi.size();

    this->Tc.resize(size);
    this->Pc.resize(size);
    this->Tr.resize(size);
    this->Pr.resize(size);
    this->Xi.resize(size);
    this->XiNew.resize(size);
    this->Alpha.resize(size);
    this->Ap.resize(size);
    this->Bp.resize(size);
    
    this->Ab.resize(size);
    this->Kij.resize(size);
    for (size_t i = 0; i < size; ++i)
    {
        this->Ab[i].resize(size);
        this->Kij[i].resize(size);
    }

    this->Ki.resize(size);
    this->M.resize(size);
    this->Ai.resize(size);
    this->Bi.resize(size);
    this->Di.resize(size);
    this->Fiv.resize(size);
    this->Fil.resize(size);

    for (size_t i = 0; i < size; ++i)
    {
        this->Tc[i] = uc.TemperatureUnits.CelciusToRankine(tc[i]);
        this->Pc[i] = uc.PressureUnits.kPaToPsi(pc[i]);
    }
    this->Pr = vc.ReducedParam(this->Pressure, this->Pc);
}

std::vector<double> DewPoint::EstimateKi(double t)
{
    std::vector<double> Result(this->Yi.size());
    for (size_t i = 0; i < this->Yi.size(); ++i)
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
    for (size_t i = 0; i < xi.size(); ++i)
    {
        Result += xi[i] * tsati[i];
    }
    return Result;
}

std::vector<double> DewPoint::EstimateTSati()
{
    std::vector<double> Result(this->Yi.size());
    for (size_t i = 0; i < Result.size(); ++i)
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
    for (size_t i = 0; i < xi.size(); ++i)
    {
        xSum += xi[i];
    }
    if (xSum != 1.0)
    {
        for (size_t i = 0; i < xi.size(); ++i)
        {
            xi[i] /= xSum;
        }
    }

    ValuesConverter vc;
    std::vector<double> Tr = vc.ReducedParam(t, this->Tc);
    this->CalculateAlpha(this->M, Tr);
    this->CalculateAp(this->Alpha, Tr, this->Pr);
    this->CalculateBp(this->Pr, Tr);
    this->CalculateAb(kij, this->Ap);
    
    this->CalculateAv();
    this->CalculateBv();
    this->CalculateZv();
    this->CalculateAl(xi);
    this->CalculateBl(xi);
    double Zl = this->CalculateZl(this->Al, this->Bl);

    std::vector<double> Fiv = this->CalculateFiv(
        this->Ab, this->Yi, this->Zv, this->Bp, this->Av, this->Bv);
    std::vector<double> Fil = this->CalculateFil(
        this->Ab, xi, Zl, this->Bp, this->Al, this->Bl);

    double s = 0.0;
    std::vector<double> XiNew(m.size());
    for (size_t i = 0; i < m.size(); ++i)
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
    this->CalculateAlpha(this->M, Tr);
    this->CalculateAp(this->Alpha, Tr, this->Pr);
    this->CalculateBp(this->Pr, Tr);
    this->CalculateAb(kij, this->Ap);
    std::vector<double> ki = this->EstimateKi(t);
    std::vector<double> xi = this->CalculateXi(ki);
    
    this->CalculateAv();
    this->CalculateBv();
    this->CalculateZv();
    this->CalculateAl(xi);
    this->CalculateBl(xi);
    double Zl = this->CalculateZl(this->Al, this->Bl);
    std::vector<double> Fiv = this->CalculateFiv(
        this->Ab, this->Yi, this->Zv, this->Bp, this->Av, this->Bv);
    std::vector<double> Fil = this->CalculateFil(
        this->Ab, xi, Zl, this->Bp, this->Al, this->Bl);

    std::vector<double> XiNew(m.size());
    for (size_t i = 0; i < m.size(); ++i)
    {
        XiNew[i] = this->Yi[i] * Fiv[i] / Fil[i];
    }

    this->Xi = xi;
    this->XiNew = XiNew;
}

double DewPoint::SelectCubicEquationRoot(
    const double &z1, const double &z2, const double &z3, 
    std::function<double (std::vector<double>)> foo)
{
    std::vector<double> roots{z1, z2, z3};
    return foo(roots);
}