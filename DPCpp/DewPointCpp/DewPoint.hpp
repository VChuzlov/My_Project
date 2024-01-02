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
};
