#include <vector>
#include <functional>

std::vector<double> VietaMethod(
    const double&, const double&, const double&);
double BrentsMethod(std::function<double (double)>, 
    double, double, double = 1.0E-8, 
    size_t MaxIter = 2000);
double Min(const std::vector<double>&);
double Max(const std::vector<double>&);
