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
      XiNew: TArrOfDouble;
      function CalculateM(const af: TArrOfDouble): TArrOfDouble;
      function CalculateAlpha(const m, tr: TArrOfDouble): TArrOfDouble;
      function CalculateAp(const alpha, tr, pr: TArrOfDouble): TArrOfDouble;
      function CalculateBp(const pr, tr: TArrOfDouble): TArrOfDouble;
      function CalculateAi(const tc, pc: TArrOfDouble): TArrOfDouble;
      function CalculateBi(const tc, pc: TArrOfDouble): TArrOfDouble;
      function CalculateDi(const m, ai, alpha, tr: TArrOfDouble): TArrOfDouble;
      function CalculateAb(
        const kij: TMatrixOfDouble; const ap: TArrOfDouble): TMatrixOfDouble;
      function CalculateAv(
        const y: TArrOfDouble; const ab: TMatrixOfDouble): Double;
      function CalculateBv(const y, bp: TArrOfDouble): Double;
      function CalculateAl(
        const x: TArrOfDouble; const ab: TMatrixOfDouble): Double;
      function CalculateBl(
        const x, bp: TArrOfDouble): Double;
      function CalculateBbl(const x, bi: TArrOfDouble): Double;
      function CalculateBbv(const y, bi: TArrOfDouble): Double;
      function CalculateAalpha(
        const mf: TArrOfDouble;
        const kij: TMatrixOfDouble;
        const ai: TArrOfDouble;
        const alpha: TArrOfDouble
      ): Double;
      function CalculateD(
        const mf, m: TArrOfDouble;
        const kij: TMatrixOfDouble;
        const ai, alpha, tr: TArrOfDouble): Double;
      function SelectCubicEquationRoot(
        const z1, z2, z3: Double; f: TFoo): Double;
      function CalculateZv(
        const av, bv: Double; method: TCubicEquationMethod): Double;
      function CalculateZl(
        const al, bl: Double; method: TCubicEquationMethod): Double;
      function CalculateFiv(
        const ab: TMatrixOfDouble;
        const y: TArrOfDouble;
        const zv: Double;
        const bp: TArrOfDouble;
        const av: Double;
        const bv: Double
      ): TArrOfDouble;
      function CalculateFil(
        const ab: TMatrixOfDouble;
        const x: TArrOfDouble;
        const zl: Double;
        const bp: TArrOfDouble;
        const al: Double;
        const bl: Double
      ): TArrOfDouble;
      function ForInitialTValue(const t: Double): Double;

    public
      constructor Create(
        const pressure: double;
        const yi: TArrOfDouble;
        const tc: TArrOfDouble;
        const pc: TArrOfDouble;
        const af: TArrOfDouble;
        const volc: TArrOfDouble
      );
      function EstimateTSati(): TArrOfDouble;
      function EstimateKi(const t: Double): TArrOfDouble;
      function CalculateXi(const ki: TArrOfDouble): TArrOfDouble;
      function EstimateTFromXiAndTSati(const xi, tsati: TArrOfDouble): Double;
      function CalculateInitialValueForT(): Double;
      function CalculateKij(
        const vc: TArrOfDouble; n: Integer = 1): TMatrixOfDouble;
      procedure PreCalculation(
        const t: Double;
        const kij: TMatrixOfDouble;
        const m: TArrOfDouble);
      function InsideJob(
        const t: Double;
        const kij: TMatrixOfDouble;
        const m: TArrOfDouble;
        const xi: TArrOfDouble): Double;
      function Condition(tol: Double = 1e-6): Boolean;
      function Calculation(): Double;
  end;


implementation

{TDewPoint}

function TDewPoint.CalculateAalpha(
  const mf: TArrOfDouble;
  const kij: TMatrixOfDouble;
  const ai, alpha: TArrOfDouble): Double;
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

function TDewPoint.CalculateAb(
  const kij: TMatrixOfDouble; const ap: TArrOfDouble): TMatrixOfDouble;
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

function TDewPoint.CalculateAi(const tc, pc: TArrOfDouble): TArrOfDouble;
var
  i: Integer;
begin
  SetLength(Result, Length(tc));
  for i := 0 to High(tc) do
    Result[i] := 0.45724 * power(8.314 * tc[i], 2) / pc[i];
end;

function TDewPoint.CalculateAl(
  const x: TArrOfDouble; const ab: TMatrixOfDouble): Double;
var
  i: Integer;
  j: Integer;
begin
  Result := 0.0;
  for i := 0 to High(x) do
    for j := 0 to High(x) do
      Result := Result + x[i] * x[j] * ab[i, j]
end;

function TDewPoint.CalculateAlpha(const m, tr: TArrOfDouble): TArrOfDouble;
var
  i: Integer;
begin
  SetLength(Result, Length(m));
  for i := 0 to High(m) do
    Result[i] := Power((1 + m[i] * (1 - Power(tr[i], 0.5))), 2);
end;

function TDewPoint.CalculateAp(const alpha, tr, pr: TArrOfDouble): TArrOfDouble;
var
  i: Integer;
begin
  SetLength(Result, Length(alpha));
  for i := 0 to High(alpha) do
    Result[i] := 0.457235529 * alpha[i] * pr[i] / Power(tr[i], 2);
end;

function TDewPoint.CalculateAv(
  const y: TArrOfDouble; const ab: TMatrixOfDouble): Double;
var
  i: Integer;
  j: Integer;
begin
  Result := 0.0;
  for i := 0 to High(y) do
    for j := 0 to High(y) do
      Result := Result + y[i] * y[j] * ab[i, j];
end;

function TDewPoint.CalculateBbl(const x, bi: TArrOfDouble): Double;
var
  i: Integer;
begin
  Result := 0.0;
  for i := 0 to High(x) do
    Result := Result + x[i] * bi[i];
end;

function TDewPoint.CalculateBbv(const y, bi: TArrOfDouble): Double;
var
  i: Integer;
begin
  Result := 0.0;
  for i := 0 to High(y) do
    Result := Result + y[i] * bi[i];
end;

function TDewPoint.CalculateBi(const tc, pc: TArrOfDouble): TArrOfDouble;
var
  i: Integer;
begin
  SetLength(Result, Length(tc));
  for i := 0 to High(tc) do
    REsult[i] := 0.07780 * 8.314 * tc[i] / pc[i];
end;

function TDewPoint.CalculateBl(
  const x, bp: TArrOfDouble): Double;
var
  i: Integer;
begin
  Result := 0.0;
  for i := 0 to High(x) do
    Result := Result + x[i] * bp[i];
end;

function TDewPoint.CalculateBp(const pr, tr: TArrOfDouble): TArrOfDouble;
var
  i: Integer;
begin
  SetLength(Result, Length(tr));
  for i := 0 to High(tr) do
    Result[i] := 0.07796074 * pr[i] / tr[i];
end;

function TDewPoint.CalculateBv(const y, bp: TArrOfDouble): Double;
var
  i: Integer;
begin
  Result := 0.0;
  for i := 0 to High(y) do
    Result := Result + y[i] * bp[i];
end;

function TDewPoint.CalculateD(
  const mf, m: TArrOfDouble;
  const kij: TMatrixOfDouble;
  const ai, alpha, tr: TArrOfDouble): Double;
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

function TDewPoint.CalculateDi(const m, ai, alpha, tr: TArrOfDouble): TArrOfDouble;
var
  i: Integer;
begin
  SetLength(Result, Length(m));
  for i := 0 to High(m) do
    Result[i] := m[i] * ai[i] * alpha[i] * Power(tr[i] / alpha[i], 0.5);
end;

function TDewPoint.CalculateFil(
  const ab: TMatrixOfDouble;
  const x: TArrOfDouble;
  const zl: Double;
  const bp: TArrOfDouble;
  const al, bl: Double): TArrOfDouble;
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

function TDewPoint.CalculateFiv(
  const ab: TMatrixOfDouble;
  const y: TArrOfDouble;
  const zv: Double;
  const bp: TArrOfDouble;
  const av, bv: Double): TArrOfDouble;
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

function TDewPoint.CalculateInitialValueForT: Double;
var
  f: TObjectiveFunction;

begin
  f :=
  function(t: Double): Double
  begin
    Result := self.ForInitialTValue(t);
  end;

  Result := Bisections(
    f, 1e2, 1200.0
  );
end;

function TDewPoint.CalculateKij(
  const vc: TArrOfDouble; n: Integer): TMatrixOfDouble;
var
  VcR3: TArrOfDouble;
  Numerator: TMatrixOfDouble;
  Denominator: TMatrixOfDouble;
  i: Integer;
  j: Integer;
begin
  SetLength(Result, Length(vc));
  SetLength(VcR3, Length(vc));
  SetLength(Numerator, Length(vc));
  SetLength(Denominator, Length(vc));
  for i := 0 to High(vc) do
    begin
      SetLength(Numerator[i], Length(vc));
      SetLength(Denominator[i], Length(vc));
      SetLength(Result[i], Length(vc));

      VcR3[i] := Power(vc[i], 1 / 3);
    end;

  for i := 0 to High(vc) do
    for j := 0 to High(vc) do
    begin
      Numerator[i, j] := Power(VcR3[i] * VcR3[j], 0.5);
      Denominator[i, j] := (VcR3[i] + VcR3[j]) / 2;
      Result[i, j] := 1 - Power(Numerator[i, j] / Denominator[i, j], n);
    end;
end;

function TDewPoint.CalculateM(const af: TArrOfDouble): TArrOfDouble;
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

function TDewPoint.CalculateXi(const ki: TArrOfDouble): TArrOfDouble;
var
  i: Integer;
begin
  SetLength(Result, Length(ki));
  for i := 0 to High(ki) do
  begin
    if ki[i] <= 1e-12 then
      ki[i] := 1e-12;
    Result[i] := self.Yi[i] / ki[i];
  end;
end;

function TDewPoint.CalculateZl(
  const al, bl: Double;
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

function TDewPoint.CalculateZv(
  const av, bv: Double;
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

function TDewPoint.Calculation: Double;
var
  T: Double;
  Kij: TMatrixOfDouble;
  m: TArrOfDouble;
  i: Integer;
  foo: TObjectiveFunction;

  uc: TUnitsConverter;
  _t: Double;
begin
  T := self.CalculateInitialValueForT();
  Kij := self.CalculateKij(self.Vc);
  m := self.CalculateM(self.Af);
  self.PreCalculation(T, Kij, m);

  foo :=
  function(t: Double): Double
  begin
    Result := self.InsideJob(t, Kij, m, self.XiNew);
  end;

  uc := TUnitsConverter.Create();
  i := 0;
  while not self.Condition() do
  begin
    i := i + 1;
    Result := Brent(
      foo, 0.8 * T,  1.2 * T
    );
    _t := uc.Temperature.RankineToCelsius(Result);
  end;

end;

function TDewPoint.Condition(tol: Double): Boolean;
var
  i: Integer;
  n: Integer;
  s: Double;
begin
  n := Length(self.XiNew);
  s := 0.0;
  for i := 0 to High(self.XiNew) do
    s := s + (Power(self.Xi[i] - self.XiNew[i], 2) / self.XiNew[i]) / n;

  Result := s <= tol;
end;

constructor TDewPoint.Create(
  const pressure: Double;
  const yi: TArrOfDouble;
  const tc: TArrOfDouble;
  const pc: TArrOfDouble;
  const af: TArrOfDouble;
  const volc: TArrOfDouble
);
var
  uc: TUnitsConverter;
  vc: TValuesConverter;
  i: Integer;
begin
  uc := TUnitsConverter.Create();
  vc := TValuesConverter.Create();
  self.Pressure := uc.Pressure.kPaToPsi(pressure);
  self.Yi := yi;
  self.Af := af;
  self.Vc := volc;

  SetLength(self.Tc, Length(tc));
  SetLength(self.Pc, Length(pc));
  SetLength(self.Tr, Length(tc));
  SetLength(self.Pr, Length(pc));
  SetLength(self.Xi, Length(yi));
  SetLength(self.XiNew, Length(yi));

  for i := 0 to High(tc) do
    begin
      self.Tc[i] := uc.Temperature.CelsiusToRankine(tc[i]);
      self.Pc[i] := uc.Pressure.kPaToPsi(pc[i]);
    end;
  self.Pr := vc.ReducedParam(self.Pressure, self.Pc);
end;


function TDewPoint.EstimateKi(const t: Double): TArrOfDouble;
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

function TDewPoint.EstimateTFromXiAndTSati(
  const xi, tsati: TArrOfDouble): Double;
var
  i: Integer;
begin
  Result := 0.0;
  for i := 0 to High(xi) do
    Result := Result + xi[i] * tsati[i];
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

function TDewPoint.ForInitialTValue(const t: Double): Double;
var
  ki: TArrOfDouble;
  xi: TArrOfDouble;
  tsati: TArrOfDouble;
  t_: Double;
  begin
    ki := self.EstimateKi(t);
    xi := self.CalculateXi(ki);
    tsati := self.EstimateTSati();
    t_ := self.EstimateTFromXiAndTSati(xi, tsati);
    Result := t - t_;
  end;


function TDewPoint.InsideJob(
  const t: Double;
  const kij: TMatrixOfDouble;
  const m, xi: TArrOfDouble): Double;
var
  vc: TValuesConverter;
  Tr: TArrOfDouble;
  Alpha: TArrOfDouble;
  Ap: TArrOfDouble;
  Bp: TArrOfDouble;
  Ab: TMatrixOfDouble;
  ki: TArrOfDouble;
  Av: Double;
  Bv: Double;
  Zv: Double;
  Al: Double;
  Bl: Double;
  Zl: Double;
  Fiv: TArrOfDouble;
  Fil: TArrOfDouble;
  XiNew: TArrOfDouble;
  s: Double;
  i: Integer;
  xSum: Double;
begin
  SetLength(XiNew, Length(m));
  xSum := 0.0;
  for i := 0 to High(xi) do
    xSum := xSum + xi[i];
  if xSum <> 1.0 then
    for i := 0 to High(xi) do
      xi[i] := xi[i] / xSum;

  vc := TValuesConverter.Create();
  Tr := vc.ReducedParam(t, self.Tc);
  Alpha := self.CalculateAlpha(m, Tr);
  Ap := self.CalculateAp(Alpha, Tr, self.Pr);
  Bp := self.CalculateBp(self.Pr, Tr);
  Ab := self.CalculateAb(kij, Ap);
  Av := self.CalculateAl(self.Yi, Ab);
  Bv := self.CalculateBv(self.Yi, Bp);
  Zv := self.CalculateZv(Av, Bv, VietaMethod);
  Al := self.CalculateAl(xi, Ab);
  Bl := self.CalculateBl(xi, Bp);
  Zl := self.CalculateZl(Al, Bl, VietaMethod);
  Fiv := self.CalculateFiv(Ab, self.Yi, Zv, Bp, Av, Bv);
  Fil := self.CalculateFil(Ab, xi, Zl, Bp, Al, Bl);

  s := 0.0;
  for i := 0 to High(m) do
  begin
    XiNew[i] := self.Yi[i] * Fiv[i] / Fil[i];
    s := s + XiNew[i];
  end;

  self.Xi := xi;
  self.XiNew := XiNew;

  Result := 1 - s;
end;

procedure TDewPoint.PreCalculation(
  const t: Double;
  const kij: TMatrixOfDouble;
  const m: TArrOfDouble);
var
  vc: TValuesConverter;
  Tr: TArrOfDouble;
  Alpha: TArrOfDouble;
  Ap: TArrOfDouble;
  Bp: TArrOfDouble;
  Ab: TMatrixOfDouble;
  ki: TArrOfDouble;
  xi: TArrOfDouble;
  Av: Double;
  Bv: Double;
  Zv: Double;
  Al: Double;
  Bl: Double;
  Zl: Double;
  Fiv: TArrOfDouble;
  Fil: TArrOfDouble;
  XiNew: TArrOfDouble;
  i: Integer;

begin
  SetLength(XiNew, Length(m));
  vc := TValuesConverter.Create();
  Tr := vc.ReducedParam(t, self.Tc);
  Alpha := self.CalculateAlpha(m, Tr);
  Ap := self.CalculateAp(Alpha, Tr, self.Pr);
  Bp := self.CalculateBp(self.Pr, Tr);
  Ab := self.CalculateAb(kij, Ap);
  ki := self.EstimateKi(t);
  xi := self.CalculateXi(ki);
  Av := self.CalculateAl(self.Yi, Ab);
  Bv := self.CalculateBv(self.Yi, Bp);
  Zv := self.CalculateZv(Av, Bv, VietaMethod);
  Al := self.CalculateAl(xi, Ab);
  Bl := self.CalculateBv(xi, Bp);
  Zl := self.CalculateZl(Al, Bl, VietaMethod);
  Fiv := self.CalculateFiv(Ab, self.Yi, Zv, Bp, Av, Bv);
  Fil := self.CalculateFil(Ab, xi, Zl, Bp, Al, Bl);

  for i := 0 to High(m) do
    XiNew[i] := self.Yi[i] * Fiv[i] / Fil[i];

  self.Xi := xi;
  self.XiNew := XiNew;

end;

function TDewPoint.SelectCubicEquationRoot(
  const z1, z2, z3: Double; f: TFoo): Double;
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
