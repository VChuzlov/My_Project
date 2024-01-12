#include <vector>

class Temperature
{
public:
    double CelciusToKelvin(double);
    double CelciusToFahreheit(double);
    double CelciusToRankine(double);
    double FahrenheitToCelcius(double);
    double FahrenheitToKelvin(double);
    double FahrenheitToRankine(double);
    double KelvinToCelcius(double);
    double KelvinToFahrenheit(double);
    double KelvinToRankine(double);
    double RankineToCelcius(double);
    double RankineToKelvin(double);
    constexpr Temperature();
};


class Pressure
{
public:
    double BarTokPa(double);
    double KgCm2TokPa(double);
    double kPaToBar(double);
    double kPaToKgCm2(double);
    double kPaToPsi(double);
    double MPaTokPa(double);
    double PsiTokPa(double);
    constexpr Pressure();
};


class UnitsConverter
{
public:
    UnitsConverter();
    Temperature TemperatureUnits;
    Pressure PressureUnits;
};

class ValuesConverter
{
public:
    std::vector<double> ReducedParam(double, std::vector<double>);
};
