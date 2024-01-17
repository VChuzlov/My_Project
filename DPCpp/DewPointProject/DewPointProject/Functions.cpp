#include <vector>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <functional>
#include "Functions.hpp"

const double Pi = 3.14159265;

std::vector<double> VietaMethod(
    const double &a, const double &b, const double &c)
{
    double 
        x1 = 0.0, 
        x2 = 0.0, 
        x3 = 0.0;
    double q, r, s;
    double fi, s_;
    std::vector<double> Roots(3);

    q = (a * a - 3. * b) / 9.;
    r = (2. * a * a * a - 9. * a * b + 27. * c) / 54.;
    s = q * q * q - r * r;

    if (s > 0.0)
    {
        fi = 1. / 3. * acos(r / pow(q, (3. / 2.)));
        x1 = -2. * sqrt(q) * cos(fi) - a / 3.;
        if (x1 < 0)
        {
            x1 = 0.0;
        }

        x2 = -2. * sqrt(q) * cos(fi + 2. / 3. * Pi) - a / 3.;
        if (x2 < 0)
        {
            x2 = 0.0;
        }

        x3 = -2. * sqrt(q) * cos(fi - 2. / 3. * Pi) - a / 3.;
        if (x3 < 0)
        {
            x3 = 0.0;
        }
    }

    if (s < 0)
    {
        if (q > 0)
        {
            fi = 1. / 3. * acosh(abs(r) / pow(q, (3. / 2.)));
            x1 = -2. * copysign(1., r) * sqrt(q) * cosh(fi) - a / 3.;
            if (x1 < 0)
            {
                x1 = 0.0;
            }
        }

        if (q < 0)
        {
            s_ = abs(r) / pow(abs(q), (3. / 2.));
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


double BrentsMethod(std::function<double (double)> f, 
    double lower, double upper, double tol, 
    size_t MaxIter)
{
    double a = lower;
	double b = upper;
	double fa = f(a);	// calculated now to save function calls
	double fb = f(b);	// calculated now to save function calls
	double fs = 0;		// initialize 

    if (!(fa * fb < 0))
	{
		// std::cout << "Signs of f(lower_bound) and f(upper_bound) must be opposites" << std::endl; // throws exception if root isn't bracketed
		return -11;
	}

    if (std::abs(fa) < std::abs(b))	// if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
	{
		std::swap(a, b);
		std::swap(fa, fb);
	}

    double c = a;			// c now equals the largest magnitude of the lower and upper bounds
	double fc = fa;			// precompute function evalutation for point c by assigning it the same value as fa
	bool mflag = true;		// boolean flag used to evaluate if statement later on
	double s = 0;			// Our Root that will be returned
	double d = 0;			// Only used if mflag is unset (mflag == false)

    for (size_t iter = 1; iter < MaxIter; iter++)
	{
		// stop if converged on root or error is less than tolerance
		if (std::abs(b - a) < tol)
		{
			// std::cout << "After " << iter << " iterations the root is: " << s << std::endl;
			return s;
		} // end if
		
		if (fa != fc && fb != fc)
		{
			// use inverse quadratic interopolation
			s =	  (a * fb * fc / ((fa - fb) * (fa - fc)))
				+ (b * fa * fc / ((fb - fa) * (fb - fc)))
				+ (c * fa * fb / ((fc - fa) * (fc - fb)));
		}
		else
		{
			// secant method
			s = b - fb * (b - a) / (fb - fa);
		}

			// checks to see whether we can use the faster converging quadratic && secant methods or if we need to use bisection
		if (((s < (3 * a + b) * 0.25) || (s > b)) ||
			(mflag && (std::abs(s - b) >= (std::abs(b - c) * 0.5))) ||
			(!mflag && (std::abs(s - b) >= (std::abs(c - d) * 0.5))) ||
			(mflag && (std::abs(b - c) < tol)) ||
			(!mflag && (std::abs(c - d) < tol)))
		{
			// bisection method
			s = (a + b) * 0.5;

			mflag = true;
		}
		else
		{
			mflag = false;
		}

		fs = f(s);	// calculate fs
		d = c;		// first time d is being used (wasnt used on first iteration because mflag was set)
		c = b;		// set c equal to upper bound
		fc = fb;	// set f(c) = f(b)

		if ( fa * fs < 0)	// fa and fs have opposite signs
		{
			b = s;
			fb = fs;	// set f(b) = f(s)
		}
		else
		{
			a = s;
			fa = fs;	// set f(a) = f(s)
		}

		if (std::abs(fa) < std::abs(fb)) // if magnitude of fa is less than magnitude of fb
		{
			std::swap(a, b);		// swap a and b
			std::swap(fa, fb);	// make sure f(a) and f(b) are correct after swap
		}

	} // end for

	// std::cout<< "The solution does not converge or iterations are not sufficient" << std::endl;
    return s;
}

double Min(const std::vector<double> &x)
{
    double minX = x[0];
    for (size_t i = 0; i < x.size(); ++i)
    {
        if (x[i] != 0.0)
        {
            minX = x[i];
            break;
        }
    }
    for (size_t i = 0; i < x.size(); ++i)
    {
        if ((minX > x[i]) && (x[i] != 0))
        {
            minX = x[i];
        }
    }
    return minX;
}

double Max(const std::vector<double> &x)
{
    double maxX = x[0];
    for (size_t i = 0; i < x.size(); ++i)
    {
        if (x[i] != 0.0)
        {
            maxX = x[i];
            break;
        }
    }
    for (size_t i = 0; i < x.size(); ++i)
    {
        if (maxX < x[i])
        {
            maxX = x[i];
        }
    }
    return maxX;
}
