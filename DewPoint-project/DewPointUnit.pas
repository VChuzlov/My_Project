unit DewPointUnit;

interface

uses
  ConvertersUnit, TypesUnit, Math, FunctionsUnit;

type

  TDewPoint = class
    private
      Pressure: double;
      Yi: TArrOfDouble;
      Tc: TArrOfDouble;
      Pc: TArrOfDouble;
      Af: TArrOfDouble;
      Vc: TArrOfDouble;
      Tr: TArrOfDouble;
      Pr: TArrOfDouble;
      Xi: TArrOfDouble;
      Xi_new: TArrOfDouble;
      function CalculateM(af: TArrOfDouble): TArrOfDouble;
      function CalculateAlpha(m, tr: TArrOfDouble): TArrOfDouble;
      function CalculateAp(alpha, tr, pr: TArrOfDouble): TArrOfDouble;
      function CalculateBp(pr, tr: TArrOfDouble): TArrOfDouble;
      function CalculateAi(tc, pc: TArrOfDouble): TArrOfDouble;
      function CalculateBi(tc, pc: TArrOfDouble): TArrOfDouble;
      function CalculateDi(m, ai, alpha, tr: TArrOfDouble): TArrOfDouble;
      function CalculateAb(
        kij: TMatrixOfDouble; ap: TArrOfDouble): TMatrixOfDouble;
      function CalculateAv(y: TArrOfDouble; ab: TMatrixOfDouble): Double;
      function CalculateBv(y, bp: TArrOfDouble): Double;
      function CalculateAl(x: TArrOfDouble; ab: TMatrixOfDouble): Double;
      function CalculateBl(x, bp: TArrOfDouble): Double;
      function CalculateBbl(x, bi: TArrOfDouble): Double;
      function CalculateBbv(y, bi: TArrOfDouble): Double;
      function CalculateAalpha(
        mf: TArrOfDouble;
        kij: TMatrixOfDouble;
        ai: TArrOfDouble;
        alpha: TArrOfDouble
      ): Double;
      function CalculateD(mf, m: TArrOfDouble; kij: TMatrixOfDouble;
        ai, alpha, tr: TArrOfDouble): Double;
      function SelectCubicEquationRoot(z1, z2, z3: Double; f: TFoo): Double;
      function CalculateZv(
        av, bv: Double; method: TCubicEquationMethod): Double;
      function CalculateZl(
        al, bl: Double; method: TCubicEquationMethod): Double;
      function CalculateFiv(
        ab: TMatrixOfDouble;
        y: TArrOfDouble;
        zv: Double;
        bp: TArrOfDouble;
        av: Double;
        bv: Double
      ): TArrOfDouble;
      function CalculateFil(
        ab: TMatrixOfDouble;
        x: TArrOfDouble;
        zl: Double;
        bp: TArrOfDouble;
        al: Double;
        bl: Double
      ): TArrOfDouble;

    public
      constructor Create(
        pressure: double;
        yi: TArrOfDouble;
        tc: TArrOfDouble;
        pc: TArrOfDouble;
        af: TArrOfDouble;
        volc: TArrOfDouble
      );
      function EstimateTSati(): TArrOfDouble;
      function EstimateKi(t: Double): TArrOfDouble;
  end;


implementation

{TDewPoint}

function TDewPoint.CalculateAalpha(mf: TArrOfDouble; kij: TMatrixOfDouble; ai,
  alpha: TArrOfDouble): Double;
var
  i: Integer;
  j: Integer;
begin
  Result := 0.0;
  for i := 0 to High(mf) do
    for j := 0 to High(mf) do
      Result := (
        Result + mf[i] * mf[j] * (1 - kij[i, j])
        * Power(ai[i] * alpha[i] * ai[j] * alpha[j], 0.5));
end;

function TDewPoint.CalculateAb(kij: TMatrixOfDouble;
  ap: TArrOfDouble): TMatrixOfDouble;
var
  i: Integer;
  j: Integer;
begin
  SetLength(Result, Length(ap));
  for i := 0 to High(ap) do
    SetLength(Result[i], Length(ap));

  for i := 0 to High(ap) do
    for j := 0 to High(ap) do
      Result[i, j] := (1 - kij[i, j]) * Power(ap[i] * ap[j], 0.5);
end;

function TDewPoint.CalculateAi(tc, pc: TArrOfDouble): TArrOfDouble;
var
  i: Integer;
begin
  SetLength(Result, Length(tc));
  for i := 0 to High(tc) do
    Result[i] := 0.45724 * power(8.314 * tc[i], 2) / pc[i];
end;

function TDewPoint.CalculateAl(x: TArrOfDouble; ab: TMatrixOfDouble): Double;
var
  i: Integer;
  j: Integer;
begin
  Result := 0.0;
  for i := 0 to High(x) do
    for j := 0 to High(x) do
      Result := Result + x[i] * x[j] * ab[i, j]
end;

function TDewPoint.CalculateAlpha(m, tr: TArrOfDouble): TArrOfDouble;
var
  i: Integer;
begin
  SetLength(Result, Length(m));
  for i := 0 to High(m) do
    Result[i] := Power((1 + m[i] * (1 - Power(tr[i], 0.5))), 2);
end;

function TDewPoint.CalculateAp(alpha, tr, pr: TArrOfDouble): TArrOfDouble;
var
  i: Integer;
begin
  SetLength(Result, Length(alpha));
  for i := 0 to High(alpha) do
    Result[i] := 0.457235529 * alpha[i] * pr[i] / Power(tr[i], 2);
end;

function TDewPoint.CalculateAv(y: TArrOfDouble; ab: TMatrixOfDouble): Double;
var
  i: Integer;
  j: Integer;
begin
  Result := 0.0;
  for i := 0 to High(y) do
    for j := 0 to High(y) do
      Result := Result + y[i] * y[j] * ab[i, j];
end;

function TDewPoint.CalculateBbl(x, bi: TArrOfDouble): Double;
var
  i: Integer;
begin
  Result := 0.0;
  for i := 0 to High(x) do
    Result := Result + x[i] * bi[i];
end;

function TDewPoint.CalculateBbv(y, bi: TArrOfDouble): Double;
var
  i: Integer;
begin
  Result := 0.0;
  for i := 0 to High(y) do
    Result := Result + y[i] * bi[i];
end;

function TDewPoint.CalculateBi(tc, pc: TArrOfDouble): TArrOfDouble;
var
  i: Integer;
begin
  SetLength(Result, Length(tc));
  for i := 0 to High(tc) do
    REsult[i] := 0.07780 * 8.314 * tc[i] / pc[i];
end;

function TDewPoint.CalculateBl(x, bp: TArrOfDouble): Double;
var
  i: Integer;
begin
  Result := 0.0;
  for i := 0 to High(x) do
    Result := Result + x[i] * bp[i];
end;

function TDewPoint.CalculateBp(pr, tr: TArrOfDouble): TArrOfDouble;
var
  i: Integer;
begin
  SetLength(Result, Length(tr));
  for i := 0 to High(tr) do
    Result[i] := 0.07796074 * pr[i] / tr[i];
end;

function TDewPoint.CalculateBv(y, bp: TArrOfDouble): Double;
var
  i: Integer;
begin
  Result := 0.0;
  for i := 0 to High(y) do
    Result := Result + y[i] * bp[i];
end;

function TDewPoint.CalculateD(mf, m: TArrOfDouble; kij: TMatrixOfDouble; ai,
  alpha, tr: TArrOfDouble): Double;
var
  i: Integer;
  j: Integer;
begin
  Result := 0.0;
  for i := 0 to High(m) do
    for j := 0 to High(m) do
      Result := (
        Result + mf[i] * mf[j] * m[j] * (1 - kij[i, j])
        * Power(ai[i] * alpha[i] * alpha[j] * tr[j], 0.5)
      );
end;

function TDewPoint.CalculateDi(m, ai, alpha, tr: TArrOfDouble): TArrOfDouble;
var
  i: Integer;
begin
  SetLength(Result, Length(m));
  for i := 0 to High(m) do
    Result[i] := m[i] * ai[i] * alpha[i] * Power(tr[i] / alpha[i], 0.5);
end;

function TDewPoint.CalculateFil(ab: TMatrixOfDouble; x: TArrOfDouble;
  zl: Double; bp: TArrOfDouble; al, bl: Double): TArrOfDouble;
var
  s: Double;
  i: Integer;
  j: Integer;
begin
  SetLength(Result, Length(x));
  for i := 0 to High(x) do
  begin
    s := 0.0;

    for j := 0 to High(x) do
      s := s + ab[i, j] * x[j];

    Result[i] := exp(
      (zl - 1) * bp[i] / bl - ln(zl - bl)
      - al / (2 * Power(2, 0.5) * bl)
      * (2 * s / al - bp[i] / bl)
      * ln((zl + (1 + Power(2, 0.5)) * bl)
            / (zl - (-1 + Power(2, 0.5)) * bl))
    );
  end;

end;

function TDewPoint.CalculateFiv(ab: TMatrixOfDouble; y: TArrOfDouble;
  zv: Double; bp: TArrOfDouble; av, bv: Double): TArrOfDouble;
var
  s: Double;
  i: Integer;
  j: Integer;
begin
  SetLength(Result, Length(y));
  for i := 0 to High(y) do
  begin
    s := 0.0;

    for j := 0 to High(y) do
      s := s + ab[i, j] * y[j];

    Result[i] := exp(
      (zv - 1) * bp[i] / bv - ln(zv - bv)
      - av / (2 * Power(2, 0.5) * bv)
      * (2 * s / av - bp[i] / bv)
      * ln((zv + (1 + Power(2, 0.5)) * bv)
            / (zv - (-1 + Power(2, 0.5)) * bv))
    );
  end;
end;

function TDewPoint.CalculateM(af: TArrOfDouble): TArrOfDouble;
var
  i: Integer;
begin
  SetLength(Result, Length(af));
  for i := 0 to High(af) do
    if af[i] <= 0.49 then
      Result[i] := 0.3796 + 1.54226 * af[i] - 0.26992 * Power(af[i], 2)
    else
      Result[i] := (
        0.379642 + 1.48503 * af[i]
        - 0.1644 * Power(af[i], 2)
        + 0.016667 * Power(af[i], 3)
      );
end;

function TDewPoint.CalculateZl(al, bl: Double;
  method: TCubicEquationMethod): Double;
var
  roots: TArrOfDouble;
begin
  roots := method(
    bl - 1,
    al - 2 * bl - 3 * Power(bl, 2),
    (-al + Power(bl, 2) + bl) * bl
  );
  Result := self.SelectCubicEquationRoot(
    roots[0],
    roots[1],
    roots[2],
    FunctionsUnit.Min
  );
end;

function TDewPoint.CalculateZv(av, bv: Double;
  method: TCubicEquationMethod): Double;
var
  roots: TArrOfDouble;
begin
  roots := method(
    bv - 1,
    av - 2 * bv - 3 * Power(bv, 2),
    (-av + Power(bv, 2) + bv) * bv
  );
  Result := self.SelectCubicEquationRoot(
    roots[0],
    roots[1],
    roots[2],
    FunctionsUnit.Max
  );
end;

constructor TDewPoint.Create(
  pressure: Double;
  yi: TArrOfDouble;
  tc: TArrOfDouble;
  pc: TArrOfDouble;
  af: TArrOfDouble;
  volc: TArrOfDouble
);
var
  uc: TUnitsConverter;
  i: Integer;
begin
  uc := TUnitsConverter.Create();
  self.Pressure := uc.Pressure.kPaToPsi(pressure);
  self.Yi := yi;
  self.Af := af;
  self.Vc := volc;

  SetLength(self.Tc, Length(tc));
  SetLength(self.Pc, Length(pc));
  SetLength(self.Tr, Length(tc));
  SetLength(self.Pr, Length(pc));
  SetLength(self.Xi, Length(yi));
  SetLength(self.Xi_new, Length(yi));

  for i := 0 to High(tc) do
    begin
      self.Tc[i] := uc.Temperature.CelsiusToRankine(tc[i]);
      self.Pc[i] := uc.Pressure.kPaToPsi(pc[i]);
    end;
end;


function TDewPoint.EstimateKi(t: Double): TArrOfDouble;
var
  i: Integer;
begin
  SetLength(Result, Length(self.Yi));
  for i := 0 to High(Result) do
    Result[i] := exp(
        ln(self.Pc[i] / self.Pressure)
        + ln(10) * (7 / 3) * (1 + self.Af[i])
        * (1 - self.Tc[i] / t)
    );
end;

function TDewPoint.EstimateTSati: TArrOfDouble;

begin
  SetLength(Result, Length(self.Yi));
  for var i := 0 to High(Result) do
    Result[i] := (
      self.Tc[i] / (1 - 3 * ln(self.Pressure / self.Pc[i])
                    / (ln(10) * (7 + 7 * self.Af[i]))
                    )
    );
end;

function TDewPoint.SelectCubicEquationRoot(z1, z2, z3: Double; f: TFoo): Double;
var
  roots: TArrOfDouble;
begin
  SetLength(roots, 3);
  roots[0] := z1;
  roots[1] := z2;
  roots[2] := z3;
  Result := f(roots);
end;

end.
