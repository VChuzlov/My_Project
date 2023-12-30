unit ConvertersUnit;

interface

uses
  TypesUnit;

type
  TTemperature = class
    function CelsiusToKelvin(c: Double): Double;
    function KelvinToCelsius(k: Double): Double;
    function CelciusToFahrenheit(c: Double): Double;
    function FahrenheitToCelsius(f: Double): Double;
    function CelsiusToRankine(c: Double): Double;
    function RankineToCelsius(r: Double): Double;
    function KelvinToRankine(k: Double): Double;
    function RankineToKelvin(r: Double): Double;
    function KelvinToFahrenheit(k: Double): Double;
    function FahrenheitToKelvin(f: Double): Double;
    function FahrenheitToRankine(f: Double): Double;
  end;

  TPressure = class
    function kPaToBar(x: Double): Double;
    function BarTokPa(x: Double): Double;
    function kPaToPsi(x: Double): Double;
    function PsiTokPa(x: Double): Double;
    function kPaToKgCm2(x: Double): Double;
    function KgCm2TokPa(x: Double): Double;
    function MPaTokPa(x: Double): Double;
  end;

  TUnitsConverter = class
    Temperature: TTemperature;
    Pressure: TPressure;
    constructor Create;
  end;

  TValuesConverter = class
    function ReducedParam(param: Double; cParam: TArrOfDouble): TArrOfDouble;
  end;

implementation


{ TTemperature }

function TTemperature.CelciusToFahrenheit(c: Double): Double;
begin
  Result := c * 9 / 5 + 32;
end;

function TTemperature.CelsiusToKelvin(c: Double): Double;
begin
  Result := c + 273.15;
end;

function TTemperature.CelsiusToRankine(c: Double): Double;
begin
  Result := c * 9 / 5 + 491.67;
end;

function TTemperature.FahrenheitToCelsius(f: Double): Double;
begin
  Result := (f - 32) * 5 / 9;
end;

function TTemperature.FahrenheitToKelvin(f: Double): Double;
begin
  Result := (f + 459.7) * 5 / 9;
end;

function TTemperature.FahrenheitToRankine(f: Double): Double;
begin
  Result := f + 459.67;
end;

function TTemperature.KelvinToCelsius(k: Double): Double;
begin
  Result := k - 273.15;
end;

function TTemperature.KelvinToFahrenheit(k: Double): Double;
begin
  Result := k * 9 / 5 - 459.7;
end;

function TTemperature.KelvinToRankine(k: Double): Double;
begin
  Result := k * 9 / 5;
end;

function TTemperature.RankineToCelsius(r: Double): Double;
begin
  Result := (r - 491.67) * 5 / 9;
end;

function TTemperature.RankineToKelvin(r: Double): Double;
begin
  Result := r * 5 / 9;
end;

{ TPressure }

function TPressure.BarTokPa(x: Double): Double;
begin
  Result := x * 100;
end;

function TPressure.KgCm2TokPa(x: Double): Double;
begin
  Result := x / 0.0101972;
end;

function TPressure.kPaToBar(x: Double): Double;
begin
  Result := x / 100;
end;

function TPressure.kPaToKgCm2(x: Double): Double;
begin
  Result := x * 0.0101972;
end;

function TPressure.kPaToPsi(x: Double): Double;
begin
  Result := x * 0.1450377377;
end;

function TPressure.MPaTokPa(x: Double): Double;
begin
  Result := x * 1000;
end;

function TPressure.PsiTokPa(x: Double): Double;
begin
  Result := x * 6.8947572932;
end;

{ TUnitsConverter }

constructor TUnitsConverter.Create;
begin
  Temperature := TTemperature.Create();
  Pressure := TPressure.Create();
end;

{ TValuesConverter }

function TValuesConverter.ReducedParam(param: Double;
  cParam: TArrOfDouble): TArrOfDouble;
var
  i: Integer;
begin
  SetLength(Result, Length(cParam));
  for i := 0 to High(cParam) do
    Result[i] := param / cParam[i];
end;


end.
