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
};


class UnitsConverter
{
public:
    UnitsConverter();
    Temperature Temperature;
    Pressure Pressure;
};
