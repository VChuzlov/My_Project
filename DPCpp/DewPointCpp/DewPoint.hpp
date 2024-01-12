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
    std::vector<double> CalculateM(std::vector<double>);
    std::vector<double> CalculateAlpha(const std::vector<double>&,
        const std::vector<double>&);
    std::vector<double> CalculateAp(const std::vector<double>&,
        const std::vector<double>&, const std::vector<double>&);
    std::vector<double> CalculateBp(const std::vector<double>&,
        const std::vector<double>&);
    std::vector<double> CalculateAi(const std::vector<double>&,
        const std::vector<double>&);
    std::vector<double> CalculateBi(const std::vector<double>&,
        const std::vector<double>&);
    std::vector<double> CalculateDi(const std::vector<double>&,
        const std::vector<double>&, const std::vector<double>&, 
        const std::vector<double>&);
    std::vector<std::vector<double>> CalculateAb(
        const std::vector<std::vector<double>>&, 
        const std::vector<double>&
    );
    double CalculateAv(const std::vector<double>&,
        const std::vector<std::vector<double>>&);
    double CalculateBv(
        const std::vector<double>&, const std::vector<double>&);
    double CalculateAl(const std::vector<double>&,
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
    double SelectCubicEquationRoot(double, double, double, 
        std::function<double (std::vector<double>)>);
    double CalculateZv(double, double,
        std::function<std::vector<double>
        (double, double, double)> = VietaMethod);
    double CalculateZl(double, double,
        std::function<std::vector<double>
        (double, double, double)> = VietaMethod);
    std::vector<double> CalculateFiv(
        const std::vector<std::vector<double>>&,
        const std::vector<double>&, double, 
        const std::vector<double>&,
        double, double);
    std::vector<double> CalculateFil(
        const std::vector<std::vector<double>>&,
        const std::vector<double>&, double, 
        const std::vector<double>&,
        double, double);
    double ForInitialTValue(double);

public:
    DewPoint(
        double, const std::vector<double>&, 
        const std::vector<double>&,
        const std::vector<double>&, 
        const std::vector<double>&, 
        const std::vector<double>&);
    std::vector<double> EstimateTSati();
    std::vector<double> EstimateKi(double);
    std::vector<double> CalculateXi(std::vector<double>);
    double EstimateTFromXiAndTSati(std::vector<double>, std::vector<double>);
    double CalculateInitialValueForT();
    std::vector<std::vector<double>> CalculateKij(
        std::vector<double>, unsigned int = 1);
    void PreCalculation(double, std::vector<std::vector<double>>,
        std::vector<double>);
    double InsideJob(double, std::vector<std::vector<double>>,
        std::vector<double>, std::vector<double>);
    bool Condition(double = 1.0E-6);
    double Calculation();
};
