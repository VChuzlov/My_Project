#include <vector>
#include <functional>

std::vector<double> VietaMethod(double, double, double);
double BrentsMethod(std::function<double (double)> f, 
    double lower, double upper, double tol, 
    unsigned int MaxIter);
