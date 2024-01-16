#include <vector>
#include <functional>
#include "Functions.hpp"

class DewPoint
{
private:
    double Pressure; 
    std::vector<double> Yi;
    std::vector<double> Tc;
    std::vector<double> Pc;
    std::vector<double> Af;
    std::vector<double> Vc;
    std::vector<double> Tr;
    std::vector<double> Pr;
    std::vector<double> Xi;
    std::vector<double> XiNew;
    std::vector<double> Alpha;
    std::vector<double> Ap;
    std::vector<double> Bp;
    std::vector<std::vector<double>> Ab;
    std::vector<std::vector<double>> Kij;
    std::vector<double> Ki;
    std::vector<double> M;
    double Av;
    double Bv;
    double Zv;
    double Al;
    double Bl;
    double Zl;
    std::vector<double> Fiv;
    std::vector<double> Fil;

    void CalculateM(const std::vector<double>&);
    void CalculateAlpha(
        const std::vector<double>&, const std::vector<double>&);
    void CalculateAp(
        const std::vector<double>&, const std::vector<double>&,
        const std::vector<double>&);
    void CalculateBp(
        const std::vector<double>&, const std::vector<double>&);
    std::vector<double> CalculateAi(
        const std::vector<double>&, const std::vector<double>&);
    std::vector<double> CalculateBi(
        const std::vector<double>&, const std::vector<double>&);
    std::vector<double> CalculateDi(
        const std::vector<double>&, const std::vector<double>&,
        const std::vector<double>&, const std::vector<double>&);
    std::vector<std::vector<double>> CalculateAb(
        const std::vector<std::vector<double>>&, 
        const std::vector<double>&);
    double CalculateAv(
        const std::vector<double>&,
        const std::vector<std::vector<double>>&);
    double CalculateBv(
        const std::vector<double>&, const std::vector<double>&);
    double CalculateAl(
        const std::vector<double>&,
        const std::vector<std::vector<double>>&);
    double CalculateBl(
        const std::vector<double>&, const std::vector<double>&);
    double CalculateBbl(
        const std::vector<double>&, const std::vector<double>&);
    double CalculateBbv(
        const std::vector<double>&, const std::vector<double>&);
    double CalculateAalpha(
        const std::vector<double>&, const std::vector<std::vector<double>>&, 
        const std::vector<double>&, const std::vector<double>&);
    double CalculateD(
        const std::vector<double>&, const std::vector<double>&,
        const std::vector<std::vector<double>>&, const std::vector<double>&,
        const std::vector<double>&, const std::vector<double>&);
    double SelectCubicEquationRoot(
        const double&, const double&, const double&, 
        std::function<double (std::vector<double>)>);
    double CalculateZv(
        const double&, const double&,
        std::function<std::vector<double>
        (double, double, double)> = VietaMethod);
    double CalculateZl(
        const double&, const double&,
        std::function<std::vector<double>
        (double, double, double)> = VietaMethod);
    std::vector<double> CalculateFiv(
        const std::vector<std::vector<double>>&,
        const std::vector<double>&, const double&, 
        const std::vector<double>&,
        const double&, const double&);
    std::vector<double> CalculateFil(
        const std::vector<std::vector<double>>&,
        const std::vector<double>&, const double&, 
        const std::vector<double>&,
        const double&, const double&);
    double ForInitialTValue(double);

public:
    DewPoint(
        const double&, const std::vector<double>&, 
        const std::vector<double>&,
        const std::vector<double>&, 
        const std::vector<double>&, 
        const std::vector<double>&);
    std::vector<double> EstimateTSati();
    std::vector<double> EstimateKi(double);
    std::vector<double> CalculateXi(const std::vector<double>&);
    double EstimateTFromXiAndTSati(
        const std::vector<double>&, const std::vector<double>&);
    double CalculateInitialValueForT();
    void CalculateKij(const std::vector<double>&, int = 1);
    void PreCalculation(
        double, const std::vector<std::vector<double>>&,
        const std::vector<double>&);
    double InsideJob(
        double, const std::vector<std::vector<double>>&,
        const std::vector<double>&, std::vector<double>&);
    bool Condition(double = 1.0E-6);
    double Calculation();
};
