unit FunctionsUnit;

interface

uses
  TypesUnit, Math;

function Max(x: TarrOfDouble): Double;
function Min(x: TarrOfDouble): Double;
function VietaMethod(a, b, c: Double): TArrOfDouble;

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

function VietaMethod(a, b, c: Double): TArrOfDouble;
var
  x1, x2, x3: Double;
  q, r, s: Double;
  fi, s_: Double;
begin
  SetLength(Result, 3);

  q := (a * a - 3 * b) / 9;
  r := (2 * a * a * a - 9 * a * b + 27 * c) / 54;
  s := q * q * q - r * r;

  if s > 0 then
  begin
    fi := 1 / 3 * ArcCos(r / Power(q, (3 / 2)));
    x1 := -2 * sqrt(q) * cos(fi) - a / 3;

    if x1 < 0 then
        x1 := 0.0;

    x2 := -2 * sqrt(q) * cos(fi + 2 / 3 * Pi) - a / 3;
    if x2 < 0 then
      x2 := 0.0;

    x3 := -2 * sqrt(q) * cos(fi - 2 / 3 * Pi) - a / 3;

    if x3 < 0 then
        x3 := 0.0;
  end;

  if s < 0 then
  begin
    if q > 0 then
    begin
      fi := 1 / 3 * ArcCosh(abs(r) / Power(q, (3 / 2)));
      x1 := -2 * sign(r) * sqrt(q) * cosh(fi) - a / 3;

      if x1 < 0 then
        x1 := 0.0;
    end;

    if q < 0 then
    begin
      s_ := abs(r) / Power(abs(q), (3 / 2));
      fi := 1 / 3 * ArcSinh(s_);
      x1 := (-2 * sign(r) * sqrt(abs(q)) * sinh(fi) - a / 3);

      if x1 < 0 then
        x1 := 0.0;
    end;

    if q = 0 then
    begin
      x1 := -Power((c - a * a * a / 27), (1 / 3)) - a / 3;

      if x1 < 0 then
        x1 := 0.0;
    end;
  end;

  if s = 0 then
  begin
    x1 := -2 * sign(r) * sqrt(q) - a / 3;

    if x1 < 0 then
        x1 := 0.0;

    x2 := sign(r) * sqrt(q) - a / 3;

    if x2 < 0 then
        x2 := 0.0;
  end;

  Result[0] := x1;
  Result[1] := x2;
  Result[2] := x3;
end;

end.
