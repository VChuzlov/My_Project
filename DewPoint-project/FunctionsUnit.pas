unit FunctionsUnit;

interface

uses
  TypesUnit;

function Max(x: TarrOfDouble): Double;
function Min(x: TarrOfDouble): Double;

implementation

function Max(x: TArrOfDouble): Double;
var
  i: Integer;
begin
  Result := x[0];
  for i := 1 to High(x) do
    if Result < x[i] then
      Result := x[i];
end;

function Min(x: TarrOfDouble): Double;
var
  i: Integer;
begin
  Result := x[0];
  for i := 1 to High(x) do
    if Result > x[i] then
      Result := x[i];
end;

end.
