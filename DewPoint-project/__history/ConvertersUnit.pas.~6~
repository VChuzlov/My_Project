unit ConvertersUnit;

interface

type
  TTemperature = class
    function CelsiusToKelvin(c: double): double;
    function KelvinToCelsius(k: double): double;
    function CelciusToFahrenheit(c: double): double;
    function FahrenheitToCelsius(f: double): double;
    function CelsiusToRankine(c: double): double;
    function RankineToCelsius(r: double): double;
    function KelvinToRankine(k: double): double;
    function RankineToKelvin(r: double): double;
    function KelvinToFahrenheit(k: double): double;
    function FahrenheitToKelvin(f: double): double;
    function FahrenheitToRankine(f: double): double;
  end;

  TPressure = class
    function kPaToBar(x: double): double;
    function BarTokPa(x: double): double;
    function kPaToPsi(x: double): double;
    function PsiTokPa(x: double): double;
    function kPaToKgCm2(x: double): double;
    function KgCm2TokPa(x: double): double;
    function MPaTokPa(x: double): double;
  end;

  TUnitsConverter = class
    Temperature: TTemperature;
    Pressure: TPressure;
    constructor Create;
  end;

implementation


{ TTemperature }

function TTemperature.CelciusToFahrenheit(c: double): double;
begin
  Result := c * 9 / 5 + 32;
end;

function TTemperature.CelsiusToKelvin(c: double): double;
begin
  Result := c + 273.15;
end;

function TTemperature.CelsiusToRankine(c: double): double;
begin
  Result := c * 9 / 5 + 491.67;
end;

function TTemperature.FahrenheitToCelsius(f: double): double;
begin
  Result := (f - 32) * 5 / 9;
end;

function TTemperature.FahrenheitToKelvin(f: double): double;
begin
  Result := (f + 459.7) * 5 / 9;
end;

function TTemperature.FahrenheitToRankine(f: double): double;
begin
  Result := f + 459.67;
end;

function TTemperature.KelvinToCelsius(k: double): double;
begin
  Result := k - 273.15;
end;

function TTemperature.KelvinToFahrenheit(k: double): double;
begin
  Result := k * 9 / 5 - 459.7;
end;

function TTemperature.KelvinToRankine(k: double): double;
begin
  Result := k * 9 / 5;
end;

function TTemperature.RankineToCelsius(r: double): double;
begin
  Result := (r - 491.67) * 5 / 9;
end;

function TTemperature.RankineToKelvin(r: double): double;
begin
  Result := r * 5 / 9;
end;

{ TPressure }

function TPressure.BarTokPa(x: double): double;
begin
  Result := x * 100;
end;

function TPressure.KgCm2TokPa(x: double): double;
begin
  Result := x / 0.0101972;
end;

function TPressure.kPaToBar(x: double): double;
begin
  Result := x / 100;
end;

function TPressure.kPaToKgCm2(x: double): double;
begin
  Result := x * 0.0101972;
end;

function TPressure.kPaToPsi(x: double): double;
begin
  Result := x * 0.1450377377;
end;

function TPressure.MPaTokPa(x: double): double;
begin
  Result := x * 1000;
end;

function TPressure.PsiTokPa(x: double): double;
begin
  Result := x * 6.8947572932;
end;

{ TUnitsConverter }

constructor TUnitsConverter.Create;
begin
  Temperature := TTemperature.Create();
  Pressure := TPressure.Create();
end;

end.
