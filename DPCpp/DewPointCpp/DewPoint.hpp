#include <vector>
#include <functional>

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
    std::vector<double> CalculateAlpha(std::vector<double>,
        std::vector<double>);
    std::vector<double> CalculateAp(std::vector<double>,
        std::vector<double>, std::vector<double>);
    std::vector<double> CalculateBp(std::vector<double>,
        std::vector<double>);
    std::vector<double> CalculateAi(std::vector<double>,
        std::vector<double>);
    std::vector<double> CalculateBi(std::vector<double>,
        std::vector<double>);
    std::vector<double> CalculateDi(std::vector<double>,
        std::vector<double>, std::vector<double>, std::vector<double>);
    std::vector<std::vector<double>> CalculateAb(
        std::vector<std::vector<double>>, std::vector<double>
    );
    double CalculateAv(std::vector<double>,
        std::vector<std::vector<double>>);
    double CalculateBv(std::vector<double>, std::vector<double>);
    double CalculateAl(std::vector<double>,
        std::vector<std::vector<double>>);
    double CalculateBl(std::vector<double>, std::vector<double>);
    double CalculateBbl(std::vector<double>, std::vector<double>);
    double CalculateBbv(std::vector<double>, std::vector<double>);
    double CalculateAalpha(std::vector<double>,
        std::vector<std::vector<double>>, std::vector<double>,
        std::vector<double>);
    double CalculateD(std::vector<double>, std::vector<double>,
        std::vector<std::vector<double>>, std::vector<double>,
        std::vector<double>, std::vector<double>);
    double SelectCubicEquationRoot(double, double, double, 
        std::function<double (std::vector<double>)>);
    double CalculateZv(double, double, std::function<double (double)>);
    double CalculateZl(double, double, std::function<double (double)>);
    std::vector<double> CalculateFiv(std::vector<std::vector<double>>,
        std::vector<double>, double, std::vector<double>,
        double, double);
    std::vector<double> CalculateFil(std::vector<std::vector<double>>,
        std::vector<double>, double, std::vector<double>,
        double, double);
    double ForinitialTValue(double);

public:
    DewPoint(double, std::vector<double>, std::vector<double>,
        std::vector<double>, std::vector<double>, std::vector<double>);
    std::vector<double> EstimateTSati();
    std::vector<double> EstimateKi(double);
    std::vector<double> EstimateXi(std::vector<double>);
    double EstimateTFromXiAndTSati(std::vector<double>, std::vector<double>);
    double CalculateInitialValueForT();
    std::vector<std::vector<double>> CalculateKij(
        std::vector<std::vector<double>>, int = 1);
};
