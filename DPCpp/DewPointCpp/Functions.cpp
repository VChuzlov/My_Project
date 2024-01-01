#include <vector>
#include<math.h>
#include "Functions.h"

#define Pi 3.14159265

std::vector<double> VietaMethod(double a, double b, double c)
{
    double x1, x2, x3;
    double q, r, s;
    double fi, s_;
    std::vector<double> Roots(3);

    q = (a * a - 3.0 * b) / 9.0;
    r = (2.0 * a * a * a - 9.0 * a * b + 27.0 * c) / 54;
    s = q * q * q - r * r;

    if (s > 0.0)
    {
        fi = 1.0 / 3.0 * acos(r / pow(q, (3 / 2)));
        x1 = -2.0 * sqrt(q) * cos(fi) - a / 3;
        if (x1 < 0)
        {
            x1 = 0.0;
        }

        x2 = -2.0 * sqrt(q) * cos(fi + 2.0 / 3.0 * Pi) - a / 3.0;
        if (x2 < 0)
        {
            x2 = 0.0;
        }

        x3 = -2.0 * sqrt(q) * cos(fi - 2.0 / 3.0 * Pi) - a / 3.0;
        if (x3 < 0)
        {
            x3 = 0.0;
        }
    }

    if (s < 0)
    {
        if (q > 0)
        {
            fi = 1.0 / 3.0 * acosh(abs(r) / pow(q, (3.0 / 2.0)));
            x1 = -2.0 * copysign(1.0, r) * sqrt(q) * cosh(fi) - a / 3.0;
            if (x1 < 0)
            {
                x1 = 0.0;
            }
        }

        if (q < 0)
        {
            s_ = abs(r) / pow(abs(q), (3.0 / 2.));
            fi = 1. / 3. * asinh(s_);
            x1 = (-2. * copysign(1., r) * sqrt(abs(q)) * sinh(fi) - a / 3.);
            if (x1 < 0)
            {
                x1 = 0.0;
            }
        }

        if (q == 0)
        {
            x1 = -pow((c - a * a * a / 27.), (1. / 3.)) - a / 3.;
            if (x1 < 0)
            {
                x1 = 0.0;
            }
        }
    }

    if (s == 0)
    {
        x1 = -2. * copysign(1., r) * sqrt(q) - a / 3.;
        if (x1 < 0)
        {
            x1 = 0.0;
        }

        x2 = copysign(1., r) * sqrt(q) - a / 3.;
        if (x2 < 0)
        {
            x2 = 0.0;
        }
    }

    Roots[0] = x1;
    Roots[1] = x2;
    Roots[2] = x3;

    return Roots;
};

