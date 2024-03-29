#include <vector>
#include <functional>

std::vector<double> VietaMethod(double, double, double);
double BrentsMethod(std::function<double (double)>, 
    double, double, double = 1.0E-8, 
    unsigned int MaxIter = 2000);
double Min(std::vector<double>);
double Max(std::vector<double>);
