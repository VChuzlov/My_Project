#include <vector>

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
    
};
