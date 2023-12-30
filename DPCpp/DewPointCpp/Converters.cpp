#include <iostream>


class Temperature
{
public:
	double CelciusToKelvin(double c) 
	{
		return c + 273.15;
	}
	double CelciusToFahreheit(double c)
	{
		return c * 9. / 5. + 32.;
	}
	double CelsiusToRankine(double c)
	{
		return c * 9. / 5. + 491.67;
	}
	double FahrenheitToCelsius(double f)
	{
		return (f - 32.) * 5. / 9.;
	}
	double FahrenheitToKelvin(double f)
	{
		return (f + 459.7) * 5. / 9.;
	}
	double FahrenheitToRankine(double f)
	{
		return f + 459.67;
	}
	double KelvinToCelsius(double k)
	{
		return k - 273.15;
	}
	double KelvinToFahrenheit(double k)
	{
		return k * 9. / 5. - 459.7;
	}
	double KelvinToRankine(double k)
	{
		return k * 9. / 5.;
	}
	double RankineToCelsius(double r)
	{
		return (r - 491.67) * 5. / 9.;
	}
	double RankineToKelvin(double r)
	{
		return r * 5. / 9.;
	}
};


class Pressure 
{
public:
	double BarTokPa(double x)
	{
		return x * 100;
	}
	double KgCm2TokPa(double x)
	{
		return x / 0.0101972;
	}
	double kPaToBar(double x)
	{
		return x / 100;
	}
	double kPaToKgCm2(double x)
	{
		return x * 0.0101972;
	}
	double kPaToPsi(double x)
	{
		return x * 0.1450377377;
	}
	double MPaTokPa(double x)
	{
		return x * 1000;
	}
	double PsiTokPa(double x)
	{
		return x * 6.8947572932;
	}
};


class UnitsConverter
{
private:
	// детали имплементации
public:
	Temperature Temperature;
	Pressure Pressure;
};