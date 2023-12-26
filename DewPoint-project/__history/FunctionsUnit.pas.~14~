unit FunctionsUnit;

interface

uses
  TypesUnit, Math;

function Max(x: TarrOfDouble): Double;
function Min(x: TarrOfDouble): Double;
function VietaMethod(a, b, c: Double): TArrOfDouble;
function Bisections(foo: TObjectiveFunction; a, b: Double;
  eps: Double = 1e-8): Double;
function Golden(foo: TObjectiveFunction; a, b: Double;
  eps: Double = 1e-8; maxit: Integer = 5000): Double;

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

function Bisections(foo: TObjectiveFunction; a, b: Double;
  eps: Double = 1e-8): Double;
var
  x: Double;

begin
  if foo(a) * foo(b) > 0 then
  begin
    Result := Infinity;
    exit
  end;

  x := (a + b) / 2;
  while abs(a - b) <= eps do
  begin
    if foo(a) * foo(x) < 0 then
      b := x
    else
      a := x;
    x := (a + b) / 2;
  end;
  Result := x;
end;

function Golden(foo: TObjectiveFunction; a, b: Double;
  eps: Double = 1e-8; maxit: Integer = 5000): Double;
var
  x1, x2: Double;
  Phi: Double;
  t: Double;
  i: Integer;

begin
  Phi := (1 + sqrt(5)) / 2;
  i := 0;

  while abs(b - a) >= eps do
  begin
    i := i + 1;
    if i > maxit then
      break;

    t := (b - a) / Phi;
    x1 := b - t;
    x2 := a + t;
    if foo(x1) >= foo(x2) then
      a := x1
    else
      b := x2;
  end;

  Result := (a + b) / 2;
end;

end.
