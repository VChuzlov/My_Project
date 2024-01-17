#include <math.h>
#include <functional>
#include "DewPoint.hpp"
//#include "Functions.hpp" циклический импорт уже есть в DewPoint.hpp
#include "Converters.hpp"
#include <iostream>
#include <omp.h>

void DewPoint::CalculateAalpha(const std::vector<double> &mf)
{
    this->Aalpha = 0.0;
    for (int i = 0; i < mf.size(); ++i)
    {
        for (int j = 0; j < mf.size(); ++j)
        {
            this->Aalpha += (
                mf[i] * mf[j] * (1 - this->Kij[i][j])
                * pow(this->Ai[i] * this->Alpha[i] 
                        * this->Ai[j] * this->Alpha[j], 0.5));
        }
    }
}

void DewPoint::CalculateAb()
{
    int size = this->Ap.size();
    #pragma omp parallel for
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            this->Ab[i][j] = (
                1 - this->Kij[i][j]) * pow(this->Ap[i] * this->Ap[j], 0.5);
        }
    }
}

void DewPoint::CalculateAi(const std::vector<double> &tc,
    const std::vector<double> &pc)
{
    #pragma omp parallel for
    for (int i = 0; i < tc.size(); ++i)
    {
        this->Ai[i] = 0.45724 * pow(8.314 * tc[i], 2) / pc[i];
    }
}

void DewPoint::CalculateAl(const std::vector<double> &x)
{
    this->Al = 0.0;
    #pragma omp parallel for
    for (int i = 0; i < x.size(); ++i)
    {
        for (int j = 0; j < x.size(); ++j)
        {
            this->Al += x[i] * x[j] * this->Ab[i][j];
        }
    }
}

void DewPoint::CalculateAlpha()
{
    #pragma omp parallel for
    for (int i = 0; i < this->M.size(); ++i)
    {
        this->Alpha[i] = pow(
            (1 + this->M[i] * (1 - pow(this->Tr[i], 0.5))), 2);
    }
}

void DewPoint::CalculateAp()
{
    #pragma omp parallel for
    for (int i = 0; i < this->Tr.size(); ++i)
    {
        this->Ap[i] = (
            0.457235529 * this->Alpha[i] * this->Pr[i] 
                / pow(this->Tr[i], 2));
    }
}

void DewPoint::CalculateAv()
{
    this->Av = 0.0;
    #pragma omp parallel for
    for (int i = 0; i < this->Yi.size(); ++i)
    {
        for (int j = 0; j < this->Yi.size(); ++j)
        {
            this->Av += this->Yi[i] * this->Yi[j] * this->Ab[i][j];
        }
    }
}

void DewPoint::CalculateBbl(
    const std::vector<double> &x)
{
    this->Bbl = 0.0;
    #pragma omp parallel for
    for (int i = 0; i < x.size(); ++i)
    {
        this->Bbl += x[i] * this->Bi[i];
    }
}

void DewPoint::CalculateBbv()
{
    this->Bbv = 0.0;
    #pragma omp parallel for
    for (int i = 0; i < this->Yi.size(); ++i)
    {
        this->Bbv += this->Yi[i] * this->Bi[i];
    }
}

void DewPoint::CalculateBi(const std::vector<double> &tc,
    const std::vector<double> &pc)
{
    #pragma omp parallel for
    for (int i = 0; i < tc.size(); ++i)
    {
        this->Bi[i] = 0.07780 * 8.314 * tc[i] / pc[i];
    }
}

void DewPoint::CalculateBl(const std::vector<double> &x)
{
    this->Bl = 0.0;
    #pragma omp parallel for
    for (int i = 0; i < x.size(); ++i)
    {
        this->Bl += x[i] * this->Bp[i];
    }
}

void DewPoint::CalculateBp()
{
    #pragma omp parallel for
    for (int i = 0; i < this->Pr.size(); ++i)
    {
        this->Bp[i] = 0.07796074 * this->Pr[i] / this->Tr[i];
    }
}

void DewPoint::CalculateBv()
{
    this->Bv = 0;
    #pragma omp parallel for
    for (int i = 0; i < this->Yi.size(); ++i)
    {
        this->Bv += this->Yi[i] * this->Bp[i];
    }
}

void DewPoint::CalculateD(
    const std::vector<double> &mf, const std::vector<double> &tr)
{
    this->D = 0.0;
    #pragma omp parallel for
    for (int i = 0; i < mf.size(); ++i)
    {
        for (int j = 0; j < mf.size(); ++j)
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
    #pragma omp parallel for
    for (int i = 0; i < tr.size(); ++i)
    {
        this->Di[i] = (
            this->M[i] * this->Ai[i] * this->Alpha[i] 
            * pow(tr[i] / this->Alpha[i], 0.5)
        );
    }
}

void DewPoint::CalculateFil(const std::vector<double> &x)
{
    double s;
    #pragma omp parallel for
    for (int i = 0; i < x.size(); ++i)
    {
        s = 0.0;
        for (int j = 0; j < x.size(); ++j)
        {
            s += this->Ab[i][j] * x[j];
        }
        this->Fil[i] = exp(
            (this->Zl - 1) * this->Bp[i] / this->Bl - log(this->Zl - this->Bl)
            - this->Al / (2 * pow(2, 0.5) * this->Bl)
            * (2 * s / this->Al - this->Bp[i] / this->Bl)
            * log((this->Zl + (1 + pow(2, 0.5)) * this->Bl)
                    / (this->Zl - (-1 + pow(2, 0.5)) * this->Bl))
        );
    }
}

void DewPoint::CalculateFiv()
{
    double s;
    #pragma omp parallel for
    for (int i = 0; i < this->Yi.size(); ++i)
    {
        s = 0.0;
        for (int j = 0; j < this->Yi.size(); ++j)
        {
            s += this->Ab[i][j] * this->Yi[j];
        }
        this->Fiv[i] = exp(
            (this->Zv - 1) * this->Bp[i] / this->Bv - log(this->Zv - this->Bv)
            - this->Av / (2 * pow(2, 0.5) * this->Bv)
            * (2 * s / this->Av - this->Bp[i] / this->Bv)
            * log((this->Zv + (1 + pow(2, 0.5)) * this->Bv)
                    / (this->Zv - (-1 + pow(2, 0.5)) * this->Bv))
        );
    }
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
    int size = vc.size();
    std::vector<double> VcR3(size);
    
    #pragma omp parallel for
    for (int i = 0; i < size; ++i)
    {
        VcR3[i] = pow(vc[i], 1. / 3.);
    }
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
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
    #pragma omp parallel for
    for (int i = 0; i < af.size(); ++i)
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

void DewPoint::CalculateXi()
{
    double k = 1e-12;
    #pragma omp parallel for
    for (int i = 0; i < this->Ki.size(); ++i)
    {
        this->Ki[i] < 1e-12 ? k = k : k = this->Ki[i];
        this->Xi[i] = this->Yi[i] / k;
    }
}

void DewPoint::CalculateZl( 
    std::function<std::vector<double> (double, double, double)> method)
{
    this->Zl = 0.0;
    std::vector<double> roots(3);
    roots = method(
        this->Bl - 1.,
        this->Al - 2. * this->Bl - 3. * pow(this->Bl, 2.),
        (-this->Al + pow(this->Bl, 2.) + this->Bl) * this->Bl
    );
    this->Zl = this->SelectCubicEquationRoot(roots[0], roots[1], roots[2], Min);
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
    const int size = this->Vc.size();
    
    UnitsConverter uc;
    int i = 0;
    double _t;

    T = this->CalculateInitialValueForT();
    this->CalculateKij(this->Vc);
    this->CalculateM(this->Af);
    this->PreCalculation(T, this->Kij, this->M);

    auto foo = [this](double t)
    {
        return this->InsideJob(t);
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
    int n = this->XiNew.size();
    double s = 0.0;
    
    #pragma omp parallel for
    for (int i = 0; i < n; ++i)
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

    int size = yi.size();

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
    for (int i = 0; i < size; ++i)
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

    #pragma omp parallel for
    for (int i = 0; i < size; ++i)
    {
        this->Tc[i] = uc.TemperatureUnits.CelciusToRankine(tc[i]);
        this->Pc[i] = uc.PressureUnits.kPaToPsi(pc[i]);
    }
    this->Pr = vc.ReducedParam(this->Pressure, this->Pc);
}

void DewPoint::EstimateKi(double t)
{
    #pragma omp parallel for
    for (int i = 0; i < this->Yi.size(); ++i)
    {
        this->Ki[i] = exp(
            log(this->Pc[i] / this->Pressure)
            + log(10) * (7. / 3.) * (1. + this->Af[i])
            * (1. - this->Tc[i] / t)
        );
    }
}

double DewPoint::EstimateTFromXiAndTSati(
    const std::vector<double> &xi, 
    const std::vector<double> &tsati)
{
    double Result = 0.0;

    #pragma omp parallel for
    for (int i = 0; i < xi.size(); ++i)
    {
        Result += xi[i] * tsati[i];
    }
    return Result;
}

std::vector<double> DewPoint::EstimateTSati()
{
    std::vector<double> Result(this->Yi.size());

    #pragma omp parallel for
    for (int i = 0; i < Result.size(); ++i)
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
    this->EstimateKi(t);
    this->CalculateXi();
    std::vector<double> tasti = this->EstimateTSati();
    double t_ = this->EstimateTFromXiAndTSati(this->Xi, tasti);
    return t - t_;
}

double DewPoint::InsideJob(const double &t)
{
    double xSum = 0.0;

    #pragma omp parallel for
    for (int i = 0; i < this->XiNew.size(); ++i)
    {
        xSum += this->XiNew[i];
        this->Xi[i] = this->XiNew[i];
    }
    if (xSum != 1.0)
    {
        for (int i = 0; i < this->XiNew.size(); ++i)
        {
            this->XiNew[i] /= xSum;
        }
    }

    ValuesConverter vc;
    this->Tr = vc.ReducedParam(t, this->Tc);
    this->CalculateAlpha();
    this->CalculateAp();
    this->CalculateBp();
    this->CalculateAb();
    
    this->CalculateAv();
    this->CalculateBv();
    this->CalculateZv();
    this->CalculateAl(this->XiNew);
    this->CalculateBl(this->XiNew);
    this->CalculateZl();

    this->CalculateFiv();
    this->CalculateFil(this->XiNew);

    double s = 0.0;

    #pragma omp parallel for
    for (int i = 0; i < this->M.size(); ++i)
    {
        this->XiNew[i] = this->Yi[i] * this->Fiv[i] / this->Fil[i];
        s += this->XiNew[i];
    }
    
    return 1.0 - s;
}

void DewPoint::PreCalculation(
    double t, const std::vector<std::vector<double>> &kij,
    const std::vector<double> &m)
{
    ValuesConverter vc;
    this->Tr = vc.ReducedParam(t, this->Tc);
    this->CalculateAlpha();
    this->CalculateAp();
    this->CalculateBp();
    this->CalculateAb();
    this->EstimateKi(t);
    this->CalculateXi();
    
    this->CalculateAv();
    this->CalculateBv();
    this->CalculateZv();
    this->CalculateAl(this->Xi);
    this->CalculateBl(this->Xi);
    this->CalculateZl();
    this->CalculateFiv();
    this->CalculateFil(this->Xi);

    #pragma omp parallel for
    for (int i = 0; i < m.size(); ++i)
    {
        this->XiNew[i] = this->Yi[i] * this->Fiv[i] / this->Fil[i];
    }
}

double DewPoint::SelectCubicEquationRoot(
    const double &z1, const double &z2, const double &z3, 
    std::function<double (std::vector<double>)> foo)
{
    std::vector<double> roots{z1, z2, z3};
    return foo(roots);
}