#include <vector>
#include <functional>

std::vector<double> VietaMethod(double, double, double);
double BrentsMethod(std::function<double (double)> f, 
    double lower, double upper, double tol = 1.0E-8, 
    unsigned int MaxIter = 2000);
