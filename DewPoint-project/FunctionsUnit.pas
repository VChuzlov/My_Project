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
function Brent(foo: TObjectiveFunction; lower, upper: Double;
  tol: Double = 1e-8; MaxIter: Integer = 5000): Double;

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
  while abs(a - b) >= eps do
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

function Brent(foo: TObjectiveFunction; lower, upper: Double;
  tol: Double = 1e-8; MaxIter: Integer = 5000): Double;
var
  a, b: Double;
  fa, fb, fs: Double;
  tmp: Double;
  c: Double; // c now equals the largest magnitude of the lower and upper bounds
  fc: Double; // precompute function evalutation for point c by assigning
             //  it the same value as fa
  mflag: boolean; // boolean flag used to evaluate if statement later on
  s: Double;  // Our Root that will be returned
  d: Double; // Only used if mflag is unset (mflag == false)
begin
  a := lower;
  b := upper;
  fa := foo(a);
  fb := foo(b);
  fs := 0.0;

  if fa * fb > 0 then
  begin
    Result := Infinity;
    exit
  end;

  // if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
  if abs(fa) < abs(b) then
  begin
    tmp := a;
    a := b;
    b := tmp;
    tmp := fa;
    fa := fb;
    fb := tmp;
  end;

  s := a;
  fc := fa;
  mflag := True;
  s := 0.0;
  d := 0.0;

  for var i := 0 to MaxIter do
  begin
    if abs(b - a) < tol then
    begin
      Result := s;
      break
    end;

    if (fa <> fc) and (fb <> fc) then
    // use inverse quadratic interopolation
      s := (a * fb * fc / ((fa - fb) * (fa - fc)))
            + (b * fa * fc / ((fb - fa) * (fb - fc)))
				    + (c * fa * fb / ((fc - fa) * (fc - fb)))
    else
    // secant method
      s := b - fb * (b - a) / (fb - fa);

    // checks to see whether we can use the faster converging 
   //  quadratic && secant methods or if we need to use bisection
    if (	( (s < (3 * a + b) * 0.25) or (s > b) ) or
				 ( mflag and (abs(s - b) >= (abs(b - c) * 0.5)) ) or
				 ( not mflag and (abs(s - b) >= (abs(c - d) * 0.5)) ) or
				 ( mflag and (abs(b - c) < tol) ) or
				 ( not mflag and (abs(c - d) < tol))	) then
    begin
      // bisection method
			s := (a + b) * 0.5;
			mflag := true;
    end

    else
      mflag := false;

    fs := foo(s);	// calculate fs
    d := c;		// first time d is being used (wasnt used on first iteration because mflag was set)
    c := b;		// set c equal to upper bound
    fc := fb;	// set f(c) = f(b)

    if (fa * fs < 0) then	 // fa and fs have opposite signs
		begin
			b := s;
			fb := fs;	 // set f(b) = f(s)
		end
    else
		begin
			a := s;
			fa := fs;	 // set f(a) = f(s)
		end;

    if (abs(fa) < abs(fb)) then // if magnitude of fa is less than magnitude of fb
		begin
      tmp := a;  // swap a and b
      a := b;
      b := tmp;
			tmp := fa;	// make sure f(a) and f(b) are correct after swap
      fa := fb;
      fb := tmp;
		end;
  end;
end;

end.
