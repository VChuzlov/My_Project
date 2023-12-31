#include <iostream>
#include <vector>
#include "Converters.h"



double Temperature::CelciusToKelvin(double c) 
{
	return c + 273.15;
}
double Temperature::CelciusToFahreheit(double c)
{
	return c * 9. / 5. + 32.;
}
double Temperature::CelciusToRankine(double c)
{
	return c * 9. / 5. + 491.67;
}
double Temperature::FahrenheitToCelcius(double f)
{
	return (f - 32.) * 5. / 9.;
}
double Temperature::FahrenheitToKelvin(double f)
{
	return (f + 459.7) * 5. / 9.;
}
double Temperature::FahrenheitToRankine(double f)
{
	return f + 459.67;
}
double Temperature::KelvinToCelcius(double k)
{
	return k - 273.15;
}
double Temperature::KelvinToFahrenheit(double k)
{
	return k * 9. / 5. - 459.7;
}
double Temperature::KelvinToRankine(double k)
{
	return k * 9. / 5.;
}
double Temperature::RankineToCelcius(double r)
{
	return (r - 491.67) * 5. / 9.;
}
double Temperature::RankineToKelvin(double r)
{
	return r * 5. / 9.;
}


double Pressure::BarTokPa(double x)
{
	return x * 100;
}
double Pressure::KgCm2TokPa(double x)
{
	return x / 0.0101972;
}
double Pressure::kPaToBar(double x)
{
	return x / 100;
}
double Pressure::kPaToKgCm2(double x)
{
	return x * 0.0101972;
}
double Pressure::kPaToPsi(double x)
{
	return x * 0.1450377377;
}
double Pressure::MPaTokPa(double x)
{
	return x * 1000;
}
double Pressure::PsiTokPa(double x)
{
	return x * 6.8947572932;
}


UnitsConverter::UnitsConverter()
{
	Temperature TemperatureUnits;
	Pressure PressureUnits;
};


std::vector<double> ValuesConverter::ReducedParam(
	double Param, std::vector<double> cParams)
{
	std::vector<double> rParams(cParams.size());

	for (int i = 0; i < cParams.size(); ++i)
	{
		rParams[i] = Param / cParams[i];
	};
	return rParams;
};
