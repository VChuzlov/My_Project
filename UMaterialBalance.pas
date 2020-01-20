unit UMaterialBalance;

interface
uses
  UPhase_Equilibrium, Dialogs, System.SysUtils, Math;

const
  NTrays = 12;
  FeedTray = 5;

type
  arrTrays = array[1..NTrays] of double;
  TArrOfDouble = array of double;
  TArrOfArrOfDouble = array of array of double;
  Tfu = function(temp, pres, Tc, Pc, om: double): double of object;

  TMatBalance = Class

    procedure WilsonCorrelation(CritT, CritP, omega: arrComp; NTrays: integer;
      Tj, Pj: arrTrays; var Kij: TArrOfArrOfDouble);
    function Wilson(Tj: double; Pj: double; Tc, Pc, om: double): double;
    function getTsat (temp, pres, Tc, Pc, om: double): double;
    function getTdp (temp, pres, Tc, Pc, om: double): double;
    function dihotomy(f:Tfu; a, b: double; Pj: double; Tc, Pc, om: double): double;
    function getTj_0(Fj: arrTrays; zf: TArrOfArrOfDouble): arrTrays;
    function getTrayMaterialBalanceError(Fj, Uj, Wj, Lj, Vj: arrTrays): arrTrays;
    function getTrayHeatBalanceError (Fj, Uj, Wj, Lj, Vj, ef_j: arrTrays;
      Hf_l, Hf_v, H_l, H_v: arrTrays): arrTrays;
    procedure InitGuessLV(Fj: arrTrays; WD: Double; LD: Double; L0: Double;
      var Lj0: arrTrays; var Vj0: arrTrays;
      var Uj: arrTrays; var Wj: arrTrays);
    procedure Gauss_Jordan(arr: TArrOfArrOfDouble; var x: TArrOfDouble);
    procedure CalculateLiquidMoleFractions(Fj, Lj, Vj, Uj, Wj: arrTrays; zij, Kij: TArrOfArrOfDouble;
      var xij: TArrOfArrOfDouble);
    function Normalize(xij: TArrOfArrOfDouble): TArrOfArrOfDouble;
    procedure Secant(a, b: double; Tj, Pj: arrTrays; xij, yij: TArrOfArrOfDouble; var rTj: arrTrays);
    procedure CalculateVaporMoleFractions(xij, Kij: TArrOfArrOfDouble; var yij: TArrOfArrOfDouble);
    function getIdealGasCompEnthalpy(Tj: arrTrays): TArrOfArrOfDouble;
    function getIdealGasCompHeatCapasity(Tj: arrTrays): TArrOfArrOfDouble;
    function getLiquidCompHeatCapasity(Tj: arrTrays): TArrOfArrOfDouble;
    function TbpVsPressure(Tbp: arrComp; P: double): arrComp;
    function getComp_dHvap(Tbp: arrComp): arrComp;
    function IntegralIdealGasCp(T1: Double; T2: arrComp): arrComp;
    function IntegralLiquidCp(T1: arrComp; T2: Double): arrComp;
    procedure CalculateEnthalpies(Tj, Pj: arrTrays; xij, yij: TArrOfArrOfDouble; var H_l, H_v: arrTrays);
    procedure RashfordRice(zf: arrComp; T, P: double; var e: double; var xf, yf: arrComp);
    procedure getEquilibrium(zf: TArrOfArrOfDouble; Tj, Pj: arrTrays;
      var xf, yf: TArrOfArrOfDouble; var ej: arrTrays);
    function getTrays_ej(Fj, Lj, Vj: arrTrays; zf, xij, yij: TArrOfArrOfDouble; Tj, Pj: arrTrays): arrTrays;
    procedure CalculateVaporFlowRates(Hf_l, Hf_v, H_l, H_v: arrTrays; Fj, Lj, Uj, Wj: arrTrays;
      Qj: arrTrays; ej: arrTrays; LD, VD: double; var Vj: arrTrays);
    procedure CalculateHeatDuties (Fj, ej, Lj, Vj, Uj, Wj, Hf_l, Hf_v, H_l, H_v: arrTrays;
      var Qc, Qr: double);
    procedure CalculateSideFlows_and_LiquidFlows(Fj, ej, Vj, Lj0: arrTrays; LD, L0: double;
      var Wj, Uj, Lj: arrTrays);
    function getErrorValue(n: integer; calcTj, calcLj, CalcVj: TArrOfArrOfDouble): double;

    function forRegularTrays(T, P: double; xi: TArrOfDouble): double;
    function forCondenser (T, P: double; zf: TArrOfDouble): double;
    procedure MatBalCalculation(Fl, Fv, Wl, Wv: arrTrays; T1, TN, P1, PN: double;
      var L, V: arrTrays);
  private
    { Private declarations }
  public
    const
      // ����������� �����������
      Tcc: arrComp = (-82.45, 32.28, 96.75, 134.9, 152, 187.2, 196.5, 234.7, 267);
      // ����������� ��������
      Pcc: arrComp = (4641, 4484, 4257, 3648, 3797, 3375, 3334, 3032, 2737);
      // ������������� ������
      omega: arrComp = (0.0115, 0.0986, 0.1524, 0.1848, 0.201, 0.2539, 0.2222, 0.3007, 0.3498);

     // Ideal Gas Comp Heat Capasity
     IdealGasCp_a: arrComp = (4.598,	1.292,	-1.009,	-0.332,	2.266,	-2.275,	-0.866,	-1.054,	-1.229);
     IdealGasCp_b: arrComp = (0.01245,	0.04254,	0.07315,	0.09189,	0.07913,	0.121,	0.1164,	0.139,	0.1615);
     IdealGasCp_c: arrComp = (0.00000286,	-0.00001657,	-0.00003789,	-0.00004409,	-0.00002647,	-0.00006519,
                              -0.00006163,	-0.00007449,	-0.0000872);
     IdealGasCp_d: arrComp = (-0.000000002703,	0.000000002081,	0.000000007678,	0.000000006915,	-0.000000000647,
                              0.00000001367,	0.00000001267,	0.00000001551,	0.00000001829);

     // Liquid Comp Heat Capasity
     LiquidCp_a: array [1..8] of double = (10.1273,	-15.3546,	3.2008,	19.7302,	-0.8949,	-0.01489,	0.2241,	-0.04342);
     LiquidCp_b: array [1..5] of double = (0.31446,	2.5346,	-2.0242,	-0.07055,	0.07264);
     LiquidCp_R: arrComp = (1.124,	1.8314,	2.4255,	2.8962,	2.8885,	3.313,	3.385,	3.812,	4.2665);
     LiquidCp_k: arrComp = (0,	0,	0,	-0.6884,	0,	-0.7643,	0,	0,	0);

     // Heat of Formation [kJ / kmole]
     dHf298: arrComp = (-74900,	-84738,	-103890,	-134590,	-126190,	-154590,	-146490,	-167290,	-187890);

     // Molar weight
     Mr: arrComp = (16.0429000854492,	30.0699005126953,	44.0970001220703,	58.1240005493164,
                    58.1240005493164,	72.1510009765625,	72.1510009765625,	86.1779022216797,	100.205001831055);

    // Normal Boiling Point
    tbpi: arrComp = (-161.525,	-88.5999969482422,	-42.1019958496094,	-11.7299865722656,
                     -0.501989746093727,	27.8780151367188,	36.0590148925781,	68.7300048828125,
                     98.4290100097656);

    // Ideal gas heat capasity
    Cpig: arrComp = (34.92,	52.49,	73.6,	96.65,	98.49,	118.9,	120.07,	142.6,	165.2);

    // Liquid heat capasity
    Cpliq: arrComp = (0,	74.48,	119.6,	129.7,	132.42,	164.85,	167.19,	197.66,	224.721);

    { Public declarations }
  End;

  function pow(x, n: double): double;
  function powZ(x: double; n: integer): double;

implementation

function pow(x, n: double): double;
begin
  Result := exp(n * ln(x))
end;

function powZ(x: double; n: integer): double;
var
  temp: integer;
  Res: double;
begin
  temp := n;
  Res := 0;
  if temp = 0 then
    Res := 1
  else
    if temp mod 2 = 0 then
      Res := powZ(x*x, Round(temp/2))
    else
      Res := x * powZ(x, temp-1);
  Result := Res
end;

procedure TMatBalance.WilsonCorrelation(CritT: arrComp; CritP: arrComp; omega: arrComp;
  NTrays: Integer; Tj: arrTrays; Pj: arrTrays; var Kij: TArrOfArrOfDouble);
var
  i, j, k: integer;
  Ps: TArrOfArrOfDouble;

begin
  SetLength(Ps, NComp, NTrays);
  for j := 0 to NTrays-1 do
    begin
      for i := 1 to NComp do
        begin
          Ps[i-1, j] := CritP[i] / 100 * exp(5.372697 * (1 + omega[i]) * (1 - (CritT[i] + 273.15) / Tj[j+1]));
          Kij[i-1, j] := Ps[i-1, j] / (Pj[j+1] * 10);
        end;
    end;
end;

function TMatBalance.Wilson(Tj: Double; Pj: Double; Tc: Double; Pc: Double; om: Double): double;
var
  Ps: double;
  i: integer;
begin
  Ps := Pc / 100 * exp(5.372697 * (1 + om) * (1 - (Tc + 273.15) / Tj));
  Result := Ps / (Pj * 10);
end;

function TMatBalance.getTsat(temp: Double; pres: Double; Tc: Double; Pc: Double; om: Double): double;
var
    i: integer;
    Ki: double;
  begin
    Ki := Wilson(temp, pres, Tc, Pc, om);
    Result := Ki - 1;
  end;

function TMatBalance.getTdp(temp: Double; pres: Double; Tc: Double; Pc: Double; om: Double): double;
var
    i: integer;
    Ki: double;
  begin
    Ki := Wilson(temp, pres, Tc, Pc, om);
    Result := 1 / Ki - 1;
  end;

function TMatBalance.dihotomy(f: Tfu; a: Double; b: Double; Pj: Double; Tc, Pc, om: double): double;
const
  eps = 1e-5;
var
  temp: double;
  
begin
  temp := (a + b) / 2;
  if f(a, Pj, Tc, Pc, om) * f(b, Pj, Tc, Pc, om) < 0 then
    begin
      repeat
        temp := (a + b) / 2;
        if f(a, Pj, Tc, Pc, om) * f(temp, Pj, Tc, Pc, om) < 0 then
          b := temp
        else
          a := temp;
      until (abs(a - b) <= eps) or (f(temp, Pj, Tc, Pc, om) = 0);
    result := temp;
    end
  else
    ShowMessage('Dihotomy, No roots');
end;

function TMatBalance.getTj_0(Fj: arrTrays; zf: TArrOfArrOfDouble): arrTrays;

var
  i, j: integer;
  Sum_zij_Fj: double;
  S_tbp_z_F: arrComp;
  Ti_sat: arrComp;
  S_zij_Fj: arrComp;
  S1, S2: double;
  Tave: double;
  Tmin: double;
  temp, pres, Tc, Pc, om: double;
begin
  Sum_zij_Fj := 0;
  s1 := 0;
  s2 := 0;

  for i := 1 to NComp do
    begin
      S_zij_Fj[i] := 0;
      S_tbp_z_F[i] := 0;
      Ti_sat[i] := dihotomy(getTsat, 90, 900, 0.195, Tcc[i], Pcc[i], omega[i]);
    end;

  for j := 1 to Ntrays do
    for i := 1 to NComp do
      begin
        Sum_zij_Fj := Sum_zij_Fj + Fj[j] * zf[i-1, j-1];
        S_zij_Fj[i] := S_zij_Fj[i] + Fj[j] * zf[i-1, j-1];
        S_tbp_z_F[i] := S_tbp_z_F[i] + Fj[j] * zf[i-1, j-1] * Ti_sat[i];
      end;

  for i := 1 to NComp do
    s1 := s1 + S_tbp_z_F[i];

  Tave := s1 / Sum_zij_Fj;

  for i := 1 to NComp do
    s2 := s2 + abs(Ti_sat[i] - Tave) * S_zij_Fj[i];

  Tmin := Tave - s2 / Sum_zij_Fj;

  for j := 1 to NTrays do
    Result[j] := Tmin + 2 * (j - 1) / NTrays * (Tave - Tmin);
end;

function TMatBalance.getTrayMaterialBalanceError(Fj, Uj, Wj, Lj, Vj: arrTrays): arrTrays;
var
  j: integer;
  input: arrTrays;
  output: arrTrays;

begin
  input[1] := Fj[1] + Vj[2];
  input[NTrays] := Fj[NTrays] + Lj[NTrays-1];
  output[NTrays] := Wj[NTrays] + Lj[NTrays] + Vj[NTrays];

  for j := 2 to NTrays-1 do
    input[j] := Fj[j] + Lj[j-1] + Vj[j+1];
  for j := 1 to NTrays-1 do
    output[j] := Uj[j] + Wj[j] + Lj[j] + Vj[j];

  for j := 1 to NTrays do
    Result[j] := input[j] - output[j];

end;

function TMatBalance.getTrayHeatBalanceError(Fj: arrTrays; Uj: arrTrays; Wj: arrTrays;
  Lj: arrTrays; Vj: arrTrays; ef_j: arrTrays; Hf_l: arrTrays; Hf_v: arrTrays;
  H_l: arrTrays; H_v: arrTrays): arrTrays;
var
  j: integer;
  input: arrTrays;
  output: arrTrays;
begin
  input[1] := (ef_j[1] * Hf_v[1] + (1 - ef_j[1]) * Hf_l[1]) * Fj[1] + Vj[2] * H_v[2];
  input[NTrays] := (ef_j[NTrays] * Hf_v[NTrays] + (1 - ef_j[NTrays]) * Hf_l[NTrays])
    * Fj[NTrays] + Lj[NTrays-1] * H_l[NTrays-1];
  output[NTrays] := Wj[NTrays] * H_v[NTrays] + Lj[NTrays] * H_l[NTrays] + Vj[NTrays] * H_v[NTrays];

  for j := 2 to NTrays-1 do
    input[j] := (ef_j[j] * Hf_v[j] + (1 - ef_j[j]) * Hf_l[j]) * Fj[j]
      + Lj[j-1] * H_l[j-1] + Vj[j+1] * H_v[j+1];
  for j := 1 to NTrays-1 do
    output[j] := (Uj[j] + Lj[j]) * H_l[j] + (Wj[j] +  Vj[j]) * H_v[j];

  for j := 1 to NTrays do
    Result[j] := input[j] - output[j];
end;

procedure TMatBalance.InitGuessLV(Fj: arrTrays; WD: Double; LD: Double; L0: Double;
  var Lj0: arrTrays; var Vj0: arrTrays;
  var Uj: arrTrays; var Wj: arrTrays);
var
  i, j: integer;
  dj: arrTrays;
  qj: arrTrays; // feed quality
  rD: double; // ��������� ����� // ������ �� 1
  rB: double; // ������� ����� // ������ �� 1
  D, B: double;
  s: double;
  TrayMaterialBalanceError: arrTrays;

  procedure Calculation(rB: double; var Fj, Uj, Wj, Lj0, Vj0: arrTrays);
  var
    j: integer;

  begin
    //Lj0[1] := L0;
    Lj0[1] := rD * LD;
    for j := 2 to NTrays-1 do
    Lj0[j] := Lj0[j-1] + Fj[j] - Uj[j];
    Lj0[NTrays] := 9.3; // ����� ��������

    for j := 2 to NTrays do
      Vj0[j] := WD + LD + Lj0[1] - Wj[j];

    for j := 1 to NTrays do
      begin
        dj[j] := (Uj[j] + Wj[j]) / (Vj0[j] + Lj0[j]);
        qj[j] := 1; //�� ���� ���� ��� ������� ����� ��������� � ����������� �� ���������
      end;

    Vj0[NTrays] := rB * Lj0[NTrays];
    for j := NTrays-1 downto 2 do
      Vj0[j] := ((1 - qj[j]) * Fj[j] + Vj0[j+1]) / (dj[j] + 1);
    for j := 2 to NTrays-1 do
      Lj0[j] := (Fj[j] + Vj0[j+1] + Lj0[j-1]) / (dj[j] + 1) - Vj0[j];

    s := 0;
    for j := 2 to NTrays-1 do
      s := s + (qj[j] + rD) * Fj[j] + (rD + 1) * Uj[j] - Wj[j];
    B := (rD * Fj[1] + (rD + 1) * Fj[NTrays] + s) / (rD + rB + 1){9.3};

    s := 0;
    for j := 2 to NTrays-1 do
      s := s + (rB + 1 - qj[j]) * Fj[j] + (rD + 1) * Uj[j] + Wj[j];
    D := ((rB + 1) * Fj[1] + rD * Fj[NTrays] + s) / (rD + rB + 1){4.5};
  end;

  function get_rB(a, b: double; Fj, Uj, Wj, Lj0, Vj0: arrTrays): double;
  const
    eps = 1e-5;
  var
    j: integer;
    rB: double;

    function f(rB: double): double;
    begin
      Calculation(rB, Fj, Uj, Wj, Lj0, Vj0);
      Result := getTrayMaterialBalanceError(Fj, Uj, Wj, Lj0, Vj0)[1];
    end;

  begin
    if f(a) * f(b) < 0 then
      repeat
        rB := (a + b) / 2;
        if f(a) * f(rB) < 0 then
          b := rB
        else
          a := rB
      until (f(rB) = 0) or (abs(a - b) <= eps)
    else
      ShowMessage('There is not solutions for rB!');
    Result := rB;
  end;

begin
  rD := {3} L0 / LD;
  rB := {1.93548387096774}3;
  Uj[1] := LD;
  Wj[1] := WD;
  Uj[NTrays] := Fj[FeedTray] - (LD + WD);  // �.�. ����� Fj
  Vj0[1] := 0; // ������ �������� �� ���� ������������
  Lj0[1] := LD * rD;

  Calculation(rB, Fj, Uj, Wj, Lj0, Vj0);
  rB := get_rB(1e-5, 1000, Fj, Uj, Wj, Lj0, Vj0);
  Calculation(rB, Fj, Uj, Wj, Lj0, Vj0);

end;

procedure TMatBalance.Gauss_Jordan(arr: TArrOfArrOfDouble; var x: TArrOfDouble);
var
  i, j, k: integer;
  d, buf: double;
begin
  for k := 1 to NTrays do
    begin
      if arr[k-1, k-1] = 0 then
        begin
          for j := 1 to NTrays+1 do
            begin
              buf:= arr[k-1, j-1];
              arr[k-1, j-1]:= arr[k, j-1];
              arr[k, j-1]:= buf;
            end;
          d:= arr[k-1, k-1];
        end
      else
        d:= arr[k-1, k-1];
      for i := 1 to NTrays+1 do
        if d <> 0 then
          arr[k-1, i-1]:= arr[k-1, i-1] / d;
      for i := 1 to NTrays do
        if i <> k then
          begin
            d:= arr[i-1, k-1];
            for j  := 1 to NTrays+1 do
              arr[i-1, j-1]:= arr[i-1, j-1] - arr[k-1, j-1] * d;
          end;
    end;
  for k := 1 to NTrays do
    x[k-1] := arr[k-1, NTrays];
end;

procedure TMatBalance.CalculateLiquidMoleFractions(Fj: arrTrays; Lj: arrTrays;
  Vj: arrTrays; Uj: arrTrays; Wj: arrTrays; zij, Kij: TArrOfArrOfDouble; var xij: TArrOfArrOfDouble);
var
  i, j: integer;
  a_j, l_j, u_j, b_j: arrTrays;
  arr: TArrOfArrOfDouble;
  x: TArrOfDouble;
  k: Integer;
begin
  SetLength(arr, NTrays, NTrays+1);
  SetLength(x, NTrays);
  l_j[1] := 0;
  for j := 1 to NTrays do
    l_j[j+1] := Lj[j];
  for k := 1 to NComp do
    begin
      for j := 1 to NTrays do
        begin
          a_j[j] := -(Lj[j] + Uj[j] + (Vj[j] + Wj[j]) * Kij[k-1, j-1]);
          b_j[j] := -Fj[j] * zij[k-1, j-1];
        end;
      for j := 1 to NTrays-1 do
        u_j[j] := Vj[j+1] * Kij[k-1, j];

      for i := 0 to NTrays-1 do
        for j := 0 to NTrays do
          arr[i, j] := 0;

      for i := 1 to NTrays do
        begin
          arr[i-1, i-1] := a_j[i];
          arr[i-1, NTrays]:= b_j[i];
        end;
      for i := 2 to NTrays do
        arr[i-1, i-2] := l_j[i];
      for i := 1 to NTrays-1 do
        arr[i-1, i] := u_j[i];

      Gauss_Jordan(arr, x);
      for j := 1 to NTrays do
        xij[k-1, j-1] := x[j-1];
    end;
end;

function TMatBalance.Normalize(xij: TArrOfArrOfDouble): TArrOfArrOfDouble;
var
  i, j: integer;
  s: TArrOfDouble;
begin
  SetLength(s, NTrays);
  SetLength(Result, NComp, NTrays);
  for j := 0 to NTrays-1 do
    s[j] := 0;

  for j := 0 to NTrays-1 do
    for i := 0 to NComp-1 do
      begin
        xij[i, j] := abs(xij[i, j]);
        s[j] := s[j] + xij[i, j];
      end;

  for j := 0 to NTrays-1 do
    for i := 0 to NComp-1 do
      if s[j] <> 0 then
        Result[i, j] := xij[i, j] / s[j];
end;

function TMatBalance.forRegularTrays(T, P: double; xi: TArrOfDouble): double;
  var
    i: integer;
    s: double;
    Ki: Double;
  begin
    s := 0;
    for i := 0 to NComp-1 do
      begin
        Ki := Wilson(T, P, Tcc[i+1], Pcc[i+1], omega[i+1]);
        s := s + Ki * xi[i];
      end;
    Result := 1 - s;
  end;

function TMatBalance.forCondenser (T, P: double; zf: TArrOfDouble): double;
  var
    i: integer;
    s: double;
    Ki: double;
  begin
    s := 0;
    for i := 0 to NComp-1 do
      begin
        Ki := Wilson(T, P, Tcc[i+1], Pcc[i+1], omega[i+1]);
        s := s + {zf[i] *} (Ki {- 1});
      end;
    Result := 1-s {- 1};
  end;

procedure TMatBalance.Secant(a: Double; b: Double; Tj, Pj: arrTrays; xij, yij: TArrOfArrOfDouble;
  var rTj: arrTrays);
type
  Tfoo = function(T, P: double; xi: TArrOfDouble): double of object;

const
  eps = 1e-5;

var
  xi: TArrOfDouble;
  i, j: integer;
  T, P: double;
  n: integer;
  temp: TArrOfDouble;
  zf: TArrOfDouble;
  tmp: double;

  procedure SecantIterations(f: Tfoo; var temp: TArrOfDouble);
  begin
    if f(temp[n-1], P, xi) * f(temp[n], P, xi) < 0 then
        begin
          repeat
            n := n + 1;
            SetLength(temp, n+1);
            temp[n] := temp[n-1] - f(temp[n-1], P, xi) * (temp[n-1] - temp[n-2])
            / (f(temp[n-1], P, xi) - f(temp[n-2], P, xi));
          until abs(temp[n] - temp[n-1]) <= eps;
        end
      else
        ShowMessage('Secant, No Roots!');
  end;

  procedure DihotomyIterations(f: Tfoo; a, b: double; var tmp: double);
  begin
    if f(a, P, xi) * f(b, P, xi) < 0 then
      begin
        repeat
          tmp := (a + b) / 2;
          if f(a, P, xi) * f(tmp, P, xi) < 0 then
            b := tmp
          else
            a := tmp
        until (abs(a - b) <= eps) or (f(tmp, P, xi) = 0);
        tmp := (a + b) / 2;
      end
    else
      ShowMessage('Dihotomy GetTemp, No Roots!');
  end;

begin
  SetLength(xi, NComp);
  SetLength(zf, NComp);

  for j := 0 to NTrays-1 do
    begin
      T := Tj[j+1];
      P := Pj[j+1];
      n := 1;
      SetLength(temp, n+1);
      temp[0] := a;
      temp[1] := b;
      for i := 0 to NComp-1 do
        xi[i] := xij[i, j];

      if j = 0 then
        begin
          for i := 0 to NComp-1 do
            zf[i] := yij[i, 1];
          //SecantIterations(forCondenser, temp);
          DihotomyIterations(forCondenser, 50, 900, tmp);
        end
      else
        //SecantIterations(forRegularTrays, temp);
        DihotomyIterations(forRegularTrays, 50, 900, tmp);

      rTj[j+1] := {temp[n]}tmp;
        for i := 0 to n-1 do
          temp[i] := 0;
      //n := n + 1;
    end;
end;

procedure TMatBalance.CalculateVaporMoleFractions(xij: TArrOfArrOfDouble; Kij: TArrOfArrOfDouble; 
  var yij: TArrOfArrOfDouble);
var
  i, j: integer;
begin
  for i := 0 to NComp-1 do
    for j := 0 to NTrays-1 do
      yij[i, j] := Kij[i, j] * xij[i, j];
end;

function TMatBalance.getIdealGasCompEnthalpy(Tj: arrTrays): TArrOfArrOfDouble; {kJ / kg}
const
// ��������� ��� ������� ��������� ���������� ����
  a: arrComp = (-12.98,	-1.7675,	39.4889,	30.903,	67.721,	64.25,	63.198,
                74.513,	71.41);
  b: arrComp = (2.36459,	1.1429,	0.395,	0.1533,	0.00854058,	-0.131798,
                -0.0117017,	-0.096697,	-0.0968949);
  c: arrComp = (-0.00213247,	-0.0003236,	0.00211409,	0.00263479,	0.00327699,
                0.003541,	0.0033164,	0.00347649,	0.003473);
  d: arrComp = (0.0000056618,	0.0000042431,	0.000000396486,	0.0000000727226,
                -0.00000110968,	-0.0000013332,	-0.0000011705,	-0.0000013212,	-0.0000013302);
  e: arrComp = (-0.00000000372476,	-0.00000000339316,	-0.000000000667176,
                -0.000000000727896,	0.000000000176646,	0.000000000251446,
                0.000000000199636,	0.000000000252365,	0.000000000255766);
  f: arrComp = (0.000000000000860896,	0.000000000000882096,	0.000000000000167936,
                0.000000000000236736,	-6.39926E-15,	-1.29576E-14,	-8.66485E-15,
                -1.34666E-14,	-1.37726E-14);
var
  i, j: integer;
begin
  // ����������� ������ ���� � ���������
  for j := 1 to NTrays do
    for i := 1 to NComp do
      Result[i-1, j-1] := a[i] + b[i] * Tj[j] + c[i] * sqr(Tj[j]) + d[i] * sqr(Tj[j]) * Tj[j]
        + e[i] * sqr(sqr(Tj[j])) + f[i] * sqr(sqr(Tj[j])) * Tj[j];
end;

function TMatBalance.getIdealGasCompHeatCapasity(Tj: arrTrays): TArrOfArrOfDouble; {cal / (mole K)}
var
  i, j: integer;

begin
  // ����������� ������ ���� � ���������
  for j := 1 to NTrays do
    for i := 1 to NComp do
      Result[i-1, j-1] := IdealGasCp_a[i] + IdealGasCp_b[i] * Tj[j] + IdealGasCp_c[i] * sqr(Tj[j])
        + IdealGasCp_d[i] * sqr(Tj[j]) * Tj[j];
end;

function TMatBalance.getLiquidCompHeatCapasity(Tj: arrTrays): TArrOfArrOfDouble; {cal / (mole K)}
var
  i, j: integer;
  IdealGasCompCp: TArrOfArrOfDouble;
begin
  SetLength(IdealGasCompCp, NComp, NTrays);
  IdealGasCompCp := getIdealGasCompHeatCapasity(Tj);
  for j := 1 to NTrays do
    for i := 1 to NComp do
      Result[i-1, j-1] := LiquidCp_a[1] * (LiquidCp_a[2] + LiquidCp_a[3] * LiquidCp_R[i]) * Tj[j] / (Tcc[i] + 273.15)
        + (LiquidCp_a[4] + LiquidCp_a[5] * LiquidCp_R[i]) * sqr(sqr(Tj[j] / (Tcc[i] + 273.15))) * Tj[j] / (Tcc[i] + 273.15)
        + LiquidCp_a[6] * sqr(LiquidCp_R[i]) / sqr(Tj[j] / (Tcc[i] + 273.15)) + LiquidCp_a[7] * LiquidCp_R[i]
        / (sqr(Tj[j] / (Tcc[i] + 273.15)) * Tj[j] / (Tcc[i] + 273.15))
        + LiquidCp_a[8] * sqr(sqr(Tj[j] / (Tcc[i] + 273.15))) * Tj[j] / (Tcc[i] + 273.15)
        + LiquidCp_k[i] * (LiquidCp_b[1] + LiquidCp_b[2] * sqr(Tj[j] / (Tcc[i] + 273.15)) +
        LiquidCp_b[3] * sqr(sqr(Tj[j] / (Tcc[i] + 273.15))) * Tj[j] / (Tcc[i] + 273.15))
        + sqr(LiquidCp_k[i]) * (LiquidCp_b[4] + LiquidCp_b[5] * sqr(Tj[j] / (Tcc[i] + 273.15))) + IdealGasCompCp[i-1, j-1];
end;

function TMatBalance.TbpVsPressure(Tbp: arrComp; P: Double{MPa}): arrComp;
var
  i: integer;
  Tnormboil: arrComp;
begin
  for i := 1 to NComp do
    Tnormboil[i] := tbpi[i] + 273.15;
  for i := 1 to NComp do
    Result[i] := 1 / (1 / (tbp[i] + 273.15) - 8.314 * ln(P / 0.101325)
      / (getComp_dHvap(Tnormboil)[i] * 4.1868));
end;

function TMatBalance.getComp_dHvap(Tbp: arrComp): arrComp; {cal / mole}
const
  R = 1.987; {cal / (mole K)}
var
  i: integer;
begin
  for i := 1 to NComp do
    Result[i] := 1.093 * R * (Tcc[i] + 273.15) * (Tbp[i] / (Tcc[i] + 273.15)
      * (ln(Pcc[i] / 101.325) - 1) / (0.930 - Tbp[i] / (Tcc[i] + 273.15)))
end;

function TMatBalance.IntegralIdealGasCp(T1: Double; T2: arrComp): arrComp;
var
  i: integer;
begin
  for i := 1 to NComp do
    Result[i] :=
      IdealGasCp_a[i] * (T2[i] - T1)
      + IdealGasCp_b[i] * (sqr(T2[i]) - sqr(T1)) / 2
      + IdealGasCp_c[i] * (powZ(T2[i], 3) - powZ(T1, 3)) / 3
      + IdealGasCp_d[i] * (powZ(T2[i], 4) - powZ(T1, 4)) / 4;
end;

function TMatBalance.IntegralLiquidCp(T1: arrComp; T2: Double): arrComp;
var
  i: integer;
  IntegralIdealGasCompCp: arrComp;
  k: arrComp;

  function get_k: arrComp;
  var
    i: integer;
    Cn, Dn: arrComp;
    Pc, Tc, Tb: arrComp;
  begin
    for i := 1 to NComp do
      begin
        Pc[i] := Pcc[i] / 101.325;
        Tc[i] := Tcc[i] + 273.15;
        Tb[i] := tbpi[i] + 273.15;
      end;
    for i := 1 to NComp do
      begin
        Cn[i] := 4.6773 + 1.8324 * liquidCp_R[i] - 0.03501 * sqr(liquidCp_R[i]);
        Dn[i] := 0.7751 * Cn[i] - 2.6354;
        Result[i] := (-ln(Pc[i]) - Cn[i] * (1 - Tc[i] / Tb[i]) + Dn[i] * ln(Tb[i]
          / Tc[i]) - 0.4218 * (1 / (Pc[i] * sqr(Tb[i] / Tc[i]))))
        / (1 - Tc[i] / Tb[i] - ln(Tb[i] / Tc[i]));
      end;
  end;

  begin
  for i := 1 to NComp do
    begin
      if T1[i] > 298 then
        IntegralIdealGasCompCp[i] := IntegralIdealGasCp(298, T1)[i]
      else
        IntegralIdealGasCompCp[i] := -IntegralIdealGasCp(298, T1)[i]
    end;
  k := get_k;
  for i := 1 to NComp do
    Result[i] :=        {
      liquidCp_a[1] * (T2 - T1[i])
      + (liquidCp_a[2] + liquidCp_a[3] * liquidCp_R[i]) * (sqr(T2) - sqr(T1[i])) / (2 * (Tcc[i] + 273.15))
      + (liquidCp_a[4] + liquidCp_a[5] * liquidCp_R[i]) * (powZ(T2, 6) - powZ(T1[i], 6)) / (6 * powZ(Tcc[i] + 273.15, 5))
      + -(liquidCp_a[6] * sqr(liquidCp_R[i]) * sqr(Tcc[i] + 273.15)) / (T2 - T1[i])
      + -(liquidCp_a[7] * liquidCp_R[i] * powZ(Tcc[i] + 273.15, 3)) / (2 * (sqr(T2) - sqr(T1[i])))
      + - (liquidCp_a[8] * powZ(Tcc[i] + 273.15, 5)) / (4 * (powZ(T2, 4) - powZ(T1[i], 4)))
      + liquidCp_k[i] * (liquidCp_b[1] * (T2 - T1[i]) + liquidCp_b[2] * (powZ(T2, 3)
        - powZ(T1[i], 3)) / (3 * sqr(Tcc[i] + 273.15)) + liquidCp_b[3] * (powZ(T2, 6) - powZ(T1[i], 6))
          / (6 * powZ(Tcc[i] + 273.15, 5)))
     + sqr(liquidCp_k[i]) * (liquidCp_b[4] * (T2 - T1[i]) + liquidCp_b[5] * (powZ(T2, 3) - powZ(T1[i], 3))
       / (3 * sqr(Tcc[i] + 273.15)))
     + IntegralIdealGasCompCp[i];
     }
     liquidCp_a[1] * (T2 - T1[i])
      + (liquidCp_a[2] + liquidCp_a[3] * liquidCp_R[i]) * (sqr(T2) - sqr(T1[i])) / (2 * (Tcc[i] + 273.15))
      + (liquidCp_a[4] + liquidCp_a[5] * liquidCp_R[i]) * (powZ(T2, 6) - powZ(T1[i], 6)) / (6 * powZ(Tcc[i] + 273.15, 5))
      + -(liquidCp_a[6] * sqr(liquidCp_R[i]) * sqr(Tcc[i] + 273.15)) / (T2 - T1[i])
      + -(liquidCp_a[7] * liquidCp_R[i] * powZ(Tcc[i] + 273.15, 3)) / (2 * (sqr(T2) - sqr(T1[i])))
      + - (liquidCp_a[8] * powZ(Tcc[i] + 273.15, 5)) / (4 * (powZ(T2, 4) - powZ(T1[i], 4)))
      + k[i] * (liquidCp_b[1] * (T2 - T1[i]) + liquidCp_b[2] * (powZ(T2, 3)
        - powZ(T1[i], 3)) / (3 * sqr(Tcc[i] + 273.15)) + liquidCp_b[3] * (powZ(T2, 6) - powZ(T1[i], 6))
          / (6 * powZ(Tcc[i] + 273.15, 5)))
     + sqr(k[i]) * (liquidCp_b[4] * (T2 - T1[i]) + liquidCp_b[5] * (powZ(T2, 3) - powZ(T1[i], 3))
       / (3 * sqr(Tcc[i] + 273.15)))
     + IntegralIdealGasCompCp[i];
end;

procedure TMatBalance.CalculateEnthalpies(Tj, Pj: arrTrays; xij: TArrOfArrOfDouble; yij: TArrOfArrOfDouble;
  var H_l: arrTrays; var H_v: arrTrays);
var
  i, j: integer;
  compH_v: TArrOfArrOfDouble;
  compH_l: TArrOfArrOfDouble;
  dHvap: arrComp;
  CompItnLiqCp: TArrOfArrOfDouble;
  CompIntIdGasCp: TArrOfArrOfDouble;
  LiquidCompHeatCapasity: TArrOfArrOfDouble;
  Tsat: arrComp;
  Tdp: arrComp;
  s: double;

begin
  SetLength(compH_v, NComp, NTrays);
  SetLength(compH_l, NComp, NTrays);
  SetLength(CompItnLiqCp, NComp, NTrays);
  SetLength(CompIntIdGasCp, NComp, NTrays);
  SetLength(LiquidCompHeatCapasity, NComp, NTrays);

  compH_v := getIdealGasCompEnthalpy(Tj);
  LiquidCompHeatCapasity := getLiquidCompHeatCapasity(Tj);
  {
  for i := 1 to NComp do
    Tsat[i] := Tbpi[i] + 273.15;
  }
  for j := 1 to Ntrays do
    begin
      Tsat := TbpVsPressure(Tbpi, Pj[j]);
      dHvap := getComp_dHvap(Tsat);
      for i := 1 to NComp do
        begin
          if Tsat[i] >= 298 then
            CompIntIdGasCp[i-1, j-1] := IntegralIdealGasCp(298, {Tdp}Tsat)[i];
          if Tsat[i] >= Tj[j] then
            CompItnLiqCp[i-1, j-1] := IntegralLiquidCp(Tsat, Tj[j])[i];
        end;
    end;

  for j := 1 to NTrays do
    begin
      s := 0;
      for i := 1 to NComp do
        s := s + {(dHf298[i] + (CompIntIdGasCp[i-1, j-1] - dHvap[i]
          + CompItnLiqCp[i-1, j-1]) * 4.1868)} {LiquidCompHeatCapasity[i-1, j-1] * 4.1868}Cpliq[i] * (Tsat[i] + Tj[j]) * xij[i-1, j-1];
      H_l[j] := s;
    end;

  for j := 1 to NTrays do
    begin
      s := 0;
      for i := 1 to NComp do
        s := s + {compH_v[i-1, j-1] * yij[i-1, j-1] * Mr[i]}
            (dHf298[i] + CompIntIdGasCp[i-1, j-1] * 4.1868) * yij[i-1, j-1];
      H_v[j] := s;
    end;
end;

procedure TMatBalance.RashfordRice(zf: arrComp; T: Double; P: Double; var e: Double;
  var xf: arrComp; var yf: arrComp);
const
  eps = 1e-5;
var
  i: integer;
  a, b: double;
  Ki: arrComp;

  function f(e: double): double;
  var
    i: integer;
    s: double;
  begin
    s := 0;
    for i := 1 to NComp do
      if (RoundTo(T, -2) <> 0) and (RoundTo(P, -2) <> 0) then
        begin
          Ki[i] := Wilson(T, P, Tcc[i], Pcc[i], omega[i]);
          s := s + zf[i] * (Ki[i] - 1) / (1 + e * (Ki[i] - 1))
        end;
    Result := s;
  end;

  function f2(zf: arrComp; T, P: double): integer;
  var
    s1, s2: double;
    i: integer;
  begin
    s1 := 0;
    s2 := 0;
    for i := 1 to NComp do
      begin
        Ki[i] := Wilson(T, P, Tcc[i], Pcc[i], omega[i]);
        s1 := s1 + zf[i] * Ki[i];
        s2 := s2 + zf[i] / Ki[i];
      end;
    if s1 < 1 then
      Result := 1
    else
      if s2 < 1 then
        Result := 2
      else
        Result := 3
  end;

begin
  a := 0;
  b := 1;
  for i := 1 to NComp do
    begin
      xf[i] := 0;
      yf[i] := 0;
    end;
  case f2(zf, T, P) of
    1://��������
      for i := 1 to NComp do
        xf[i] := zf[i];
    2://���
      for i := 1 to NComp do
        yf[i] := zf[i];
    3:// 2 ����
      if f(a) * f(b) < 0 then
        begin
          repeat
            e := (a + b) / 2;
            if f(a) * f(e) < 0 then
              b := e
            else
              a := e
          until (abs(a - b) <= eps) or (f(e) = 0);
          for i := 1 to NComp do
            begin
              xf[i] := zf[i] / (1 + e * (Ki[i] - 1));
              yf[i] := {Ki[i] * xf[i]} zf[i] * Ki[i] / (1 + e * (Ki[i] - 1));
            end;
        end
      else
        ShowMessage('Rashford-Rice, No Roots!');
  end;

end;

procedure TMatBalance.getEquilibrium(zf: TArrOfArrOfDouble; Tj: arrTrays; Pj: arrTrays;
  var xf: TArrOfArrOfDouble; var yf: TArrOfArrOfDouble; var ej: arrTrays);
var
  i, j, k: integer;
  tray_zf: arrComp;
  tray_xf: arrComp;
  tray_yf: arrComp;
begin
  SetLength(xf, NComp, NTrays);
  SetLength(yf, NComp, NTrays);

  for j := 1 to NTrays do
    begin
      for i := 1 to NComp do
        begin
          tray_zf[i] := zf[i-1, j-1];
        end;
      RashfordRice(tray_zf, Tj[j], Pj[j], ej[j], tray_xf, tray_yf);
      for k := 1 to NComp do
        begin
          xf[k-1, j-1] := tray_xf[k];
          yf[k-1, j-1] := tray_yf[k];
        end;
    end;
  xf := Normalize(xf);
  yf := Normalize(yf);
end;

function TMatBalance.getTrays_ej(Fj: arrTrays; Lj: arrTrays; Vj: arrTrays;
  zf: TArrOfArrOfDouble; xij: TArrOfArrOfDouble; yij: TArrOfArrOfDouble; Tj, Pj: arrTrays): arrTrays;
var
  i, j, k: integer;
  z_tr: TArrOfArrOfDouble;
  s: double;
  x_tr, y_tr: TArrOfArrOfDouble;
begin
  SetLength(z_tr, NComp, NTrays);
  SetLength(x_tr, NComp, NTrays);
  SetLength(y_tr, NComp, NTrays);

  for i := 1 to NComp do
    begin
      if (Fj[1] + Vj[2]) <> 0 then
        z_tr[i-1, 0] := (Fj[1] * zf[i-1, 0] + Vj[2] * yij[i-1, 1]) / (Fj[1] + Vj[2]);
      if (Fj[NTrays] + Lj[NTrays-1]) <> 0then
        z_tr[i-1, NTrays-1] := (Fj[NTrays] * zf[i-1, NTrays-1] + Lj[NTrays-1] * xij[i-1, NTrays-2])
          / (Fj[NTrays] + Lj[NTrays-1]);
    end;

  for j := 2 to NTrays-1 do
    begin
      s := 0;
      for i := 1 to NComp do
        if (Fj[j] + Lj[j-1] + Vj[j+1]) <> 0 then
          z_tr[i-1, j-1] := (Fj[j] * zf[i-1, j-1] + Lj[j-1] * xij[i-1, j-2] + Vj[j+1] * yij[i-1, j])
            / (Fj[j] + Lj[j-1] + Vj[j+1]);
    end;
  getEquilibrium(z_tr, Tj, Pj, x_tr, y_tr, Result);
end;

procedure TMatBalance.CalculateVaporFlowRates(Hf_l: arrTrays; Hf_v: arrTrays;
  H_l: arrTrays; H_v: arrTrays; Fj: arrTrays; Lj: arrTrays; Uj: arrTrays; Wj: arrTrays;
  Qj: arrTrays; ej: arrTrays; LD, VD: double; var Vj: arrTrays);
var
  j, k: integer;
  alp, bet, gam: arrTrays;
  s: double;

begin
  Vj[NTrays] := Fj[NTrays] + Lj[NTrays-1] - Uj[NTrays] - Wj[NTrays];
  //Vj[2] := Lj[1] + LD + VD;
  for j := NTrays-1 downto 2 do
    begin
      s := 0;
      for k := 1 to j-1 do
        s := s + (Fj[k] - Uj[k] - Wj[k]);
      gam[j] := Hf_l[j] * (1 - ej[j]) * Fj[j] + Hf_v[j] * ej[j] * Fj[j] - H_l[j] * Fj[j]
        + (H_l[j-1] - H_l[j]) * s + (H_l[j] - H_v[j]) * Wj[j] - Qj[j];
      bet[j] := H_v[j+1] - H_l[j];
      alp[j] := H_v[j] - H_l[j-1];
      Vj[j] := (bet[j] * Vj[j+1] + gam[j]) / alp[j];
    end;
end;

procedure TMatBalance.CalculateHeatDuties(Fj, ej: arrTrays; Lj: arrTrays; Vj: arrTrays;
  Uj: arrTrays; Wj: arrTrays; Hf_l: arrTrays; Hf_v: arrTrays; H_l: arrTrays;
  H_v: arrTrays; var Qc: Double; var Qr: Double);
begin
  Qc := (1 - ej[1]) *Fj[1] * Hf_l[1] +  ej[1] * Fj[1] * Hf_v[1] + Vj[2] * H_v[2] - Wj[1] * H_v[1]
    - (Lj[1] + Uj[1]) * H_l[1];
  Qr := (1 - ej[NTrays]) * Fj[NTrays] * Hf_l[NTrays] + ej[NTrays] * Fj[NTrays] * Hf_v[NTrays]
    + Lj[Ntrays-1] * H_l[NTrays-1] - (Vj[NTrays] + Wj[NTrays]) * H_v[NTrays] - Uj[NTrays] * H_l[NTrays];
end;

procedure TMatBalance.CalculateSideFlows_and_LiquidFlows(Fj: arrTrays; ej: arrTrays;
  Vj: arrTrays; Lj0: arrTrays; LD, L0: double; var Wj: arrTrays; var Uj: arrTrays; var Lj: arrTrays);
var
  dj: arrTrays;
  Rj: arrTrays;
  j: integer;
  rD, rB: double;
  B, D: double;

  procedure Calculation(rB: double; var Fj, Uj, Wj, Lj, Vj: arrTrays);
  var
    j: integer;
    qj: arrTrays;
    WD: double; // ����������
    s: double;

  begin
    //Lj0[1] := L0;
    Lj[1] := rD * LD;
    for j := 2 to NTrays-1 do
    Lj[j] := Lj[j-1] + Fj[j] - Uj[j];
    Lj[NTrays] := 9.3; // ����� ��������

    for j := 2 to NTrays do
      Vj[j] := WD + LD + Lj[1] - Wj[j];

    for j := 1 to NTrays do
      begin
        dj[j] := (Uj[j] + Wj[j]) / (Vj[j] + Lj[j]);
        qj[j] := 1; //�� ���� ���� ��� ������� ����� ��������� � ����������� �� ���������
      end;

    Vj[NTrays] := rB * Lj[NTrays];
    for j := NTrays-1 downto 2 do
      Vj[j] := ((1 - qj[j]) * Fj[j] + Vj[j+1]) / (dj[j] + 1);
    for j := 2 to NTrays-1 do
      Lj[j] := (Fj[j] + Vj[j+1] + Lj[j-1]) / (dj[j] + 1) - Vj[j];

    s := 0;
    for j := 2 to NTrays-1 do
      s := s + (qj[j] + rD) * Fj[j] + (rD + 1) * Uj[j] - Wj[j];
    B := (rD * Fj[1] + (rD + 1) * Fj[NTrays] + s) / (rD + rB + 1){9.3};

    s := 0;
    for j := 2 to NTrays-1 do
      s := s + (rB + 1 - qj[j]) * Fj[j] + (rD + 1) * Uj[j] + Wj[j];
    D := ((rB + 1) * Fj[1] + rD * Fj[NTrays] + s) / (rD + rB + 1){4.5};
  end;

  function get_rB(a, b: double; Fj, Uj, Wj, Lj, Vj: arrTrays): double;
  const
    eps = 1e-5;
  var
    j: integer;
    rB: double;

    function f(rB: double): double;
    begin
      Calculation(rB, Fj, Uj, Wj, Lj, Vj);
      Result := getTrayMaterialBalanceError(Fj, Uj, Wj, Lj, Vj)[1];
    end;

  begin
    if f(a) * f(b) < 0 then
      repeat
        rB := (a + b) / 2;
        if f(a) * f(rB) < 0 then
          b := rB
        else
          a := rB
      until (f(rB) = 0) or (abs(a - b) <= eps)
    else
      ShowMessage('There is not solutions for rB!');
    Result := rB;
  end;

begin
  for j := 1 to NTrays do
    begin
      if Lj0[j] + Vj[j] <> 0 then
        dj[j] := (Uj[j] + Wj[j]) / (Lj0[j] + Vj[j])
      else
        dj[j] := 0;
      if j = 1 then
        Rj[j] := (dj[j] / (dj[j] + 1)) * (Fj[j] + Vj[j+1] {+ Lj0[j-1]})
      else
        Rj[j] := (dj[j] / (dj[j] + 1)) * (Fj[j] + Vj[j+1] + Lj0[j-1]);
      Wj[j] := ej[j] * Rj[j];
      Uj[j] := (1 - ej[j]) * Rj[j];
    end;

  Uj[1] := 4.5;
  Uj[NTrays] := 9.3;

  for j := 2 to NTrays-1 do
    Lj[j] := Fj[j] + Vj[j+1] - Vj[j] + Lj[j-1] - Rj[j];
  Lj[NTrays] := Uj[NTrays];
  Lj[1] := L0;

  rD := Lj[1] / LD;
  rB := Vj[NTrays] / Uj[NTrays];
  Calculation(rB, Fj, Uj, Wj, Lj, Vj);
  rB := get_rB(1e-5, 1000, Fj, Uj, Wj, Lj, Vj);
  Calculation(rB, Fj, Uj, Wj, Lj, Vj);

end;

function TMatBalance.getErrorValue(n: integer; calcTj, calcLj, calcVj: TArrOfArrOfDouble): double;
var
  j: integer;
  sT, sL, sV: double;
begin
  sT := 0;
  sL := 0;
  sV := 0;
  for j := 0 to NTrays-1 do
    sT := sT + sqr((calcTj[n-1, j] - calcTj[n-2, j]) / calcTj[n-1, j]);
  for j := 0 to NTrays-2 do
    sL := sL + sqr((calcLj[n-1, j] - calcLj[n-2, j]) / calcLj[n-1, j]);
  for j := 1 to NTrays-1 do
    sV := sV + sqr((calcVj[n-1, j] - calcVj[n-2, j]) / calcVj[n-1, j]);
  Result := sT + sL + sV;
end;

procedure TMatBalance.MatBalCalculation(Fl: arrTrays; Fv: arrTrays; Wl: arrTrays;
  Wv: arrTrays; T1: Double; TN: Double; P1: Double; PN: Double; var L: arrTrays; var V: arrTrays);
const
  tolerance = 1e-5;
var
  zf: TArrOfArrOfDouble;
  Res: TArrOfDouble;
  Fj, Wj, Uj: arrTrays;
  Lj0, Vj0: arrTrays;
  Tj_0: arrTrays;
  Pj: arrTrays;
  i, j: Integer;
  Kij: TArrOfArrOfDouble;
  xij: TArrOfArrOfDouble;
  norm_xij: TArrOfArrOfDouble;
  Tj: arrTrays;
  Tfj: arrTrays;
  Pfj: arrTrays;
  yij: TArrOfArrOfDouble;
  H_l, H_v: arrTrays;
  xf, yf: TArrOfArrOfDouble;
  ej: arrTrays;
  ej_tr: arrTrays;
  Hf_l, Hf_v: arrTrays;
  Lj: arrTrays;
  Vj: arrTrays;
  Qj: arrTrays;
  Qc, Qr: Double;
  calcTj: TArrOfArrOfDouble;
  calcLj: TArrOfArrOfDouble;
  calcVj: TArrOfArrOfDouble;
  n: integer;
  omegaj: ArrTrays; // ���� ���� � ������� �������
  TrayMaterialBalanceError: arrTrays;
  TrayHeatBalanceError: arrTrays;

begin
  n := 1;

  SetLength(zf, NComp, NTrays);
  SetLength(Kij, NComp, NTrays);
  SetLength(xij, NComp, NTrays);
  SetLength(norm_xij, NComp, NTrays);
  SetLength(yij, NComp, NTrays);
  SetLength(calcTj, n, NTrays);
  SetLength(calcLj, n, NTrays);
  SetLength(calcVj, n, NTrays);

  for i := 0 to NComp-1 do
    for j := 0 to NTrays-1 do
      zf[i, j] := 0;
  zf[0, FeedTray-1] := 2.22222222222222e-002;
  zf[1, FeedTray-1] := 4.44444444444444e-002;
  zf[2, FeedTray-1] := 6.66666666666667e-002;
  zf[3, FeedTray-1] := 8.88888888888889e-002;
  zf[4, FeedTray-1] := 0.111111111111111;
  zf[5, FeedTray-1] := 0.133333333333333;
  zf[6, FeedTray-1] := 0.155555555555556;
  zf[7, FeedTray-1] := 0.177777777777778;
  zf[8, FeedTray-1] := 0.200000000000000;

  for j := 1 to NTrays do
    begin
      Fj[j] := 0;
      Wj[j] := 0;
      Uj[j] := 0;
      Qj[j] := 0;
    end;
  Fj[FeedTray] := 13.8;
  Pj[1] := P1;
  Pj[Ntrays] := PN;
  for j := 2 to NTrays-1 do
    Pj[j] := Pj[j-1] + (Pj[NTrays] - Pj[1]) / NTrays;

  for j := 1 to NTrays do
    begin
      Tfj[j] := 50 + 273.15;
      Pfj[j] := Pj[j];
      omegaj[j] := 0; // �.�. ��� ����� � ���� ��������
    end;

  Tj_0 := getTj_0(Fj, zf);
  InitGuessLV(Fj, 0, 4.5, 13.5, Lj0, Vj0, Uj, Wj);
  WilsonCorrelation(Tcc, Pcc, omega, NTrays, Tj_0, Pj, Kij);
  CalculateLiquidMoleFractions(Fj, Lj0, Vj0, Uj, Wj, zf, Kij, xij);
  xij := Normalize(xij);
  Secant(90, 400, Tj_0, Pj, xij, yij, Tj);
  WilsonCorrelation(Tcc, Pcc, omega, NTrays, Tj, Pj, Kij);
  CalculateVaporMoleFractions(xij, Kij, yij);
  yij := Normalize(yij);
  CalculateEnthalpies(Tj, Pj, xij, yij, H_l, H_v);
  getEquilibrium(zf, Tfj, Pj, xf, yf, ej);
  CalculateEnthalpies(Tfj, Pfj, xf, yf, Hf_l, Hf_v);
  CalculateHeatDuties(Fj, ej, Lj0, Vj, Uj, Wj, Hf_l, Hf_v, H_l, H_v, Qc, Qr);
  TrayMaterialBalanceError := getTrayMaterialBalanceError(Fj, Uj, Wj, Lj0, Vj0);
  TrayHeatBalanceError := getTrayHeatBalanceError(Fj, Uj, Wj, Lj0, Vj0, ej, Hf_l, Hf_v, H_l, H_v);
  for j := 0 to NTrays-1 do
    begin
      calcTj[n-1, j] := Tj_0[j+1];
      calcLj[n-1, j] := Lj0[j+1];
      calcVj[n-1, j] := Vj0[j+1];
    end;
  WilsonCorrelation(Tcc, Pcc, omega, NTrays, Tj_0, Pj, Kij);
  repeat
    CalculateLiquidMoleFractions(Fj, Lj0, Vj0, Uj, Wj, zf, Kij, xij);
    xij := Normalize(xij);
    Secant(90, 400, Tj_0, Pj, xij, yij, Tj);
    WilsonCorrelation(Tcc, Pcc, omega, NTrays, Tj, Pj, Kij);
    CalculateVaporMoleFractions(xij, Kij, yij);
    yij := Normalize(yij);
    CalculateEnthalpies(Tj, Pj, xij, yij, H_l, H_v);
    getEquilibrium(zf, Tfj, Pj, xf, yf, ej);
    CalculateEnthalpies(Tfj, Pfj, xf, yf, Hf_l, Hf_v);
    CalculateVaporFlowRates(Hf_l, Hf_v, H_l, H_v, Fj, Lj0, Uj, Wj, Qj, ej, 4.5, 0, Vj);
    CalculateHeatDuties(Fj, ej, Lj0, Vj, Uj, Wj, Hf_l, Hf_v, H_l, H_v, Qc, Qr);
    ej_tr := getTrays_ej(Fj, Lj, Vj, zf, xij, yij, Tj, Pj); // ��� ���������, ��� �� �����
    CalculateSideFlows_and_LiquidFlows(Fj, omegaj, Vj, Lj0, 4.5, 13.5, Wj, Uj, Lj);
    TrayMaterialBalanceError := getTrayMaterialBalanceError(Fj, Uj, Wj, Lj, Vj);
    TrayHeatBalanceError := getTrayHeatBalanceError(Fj, Uj, Wj, Lj, Vj, ej, Hf_l, Hf_v, H_l, H_v);
    n := n + 1;
    if n >= 1e5 then
      begin
        ShowMessage('Main Calculation, No Solutions!');
        break
      end;
    SetLength(calcTj, n, NTrays);
    SetLength(calcLj, n, NTrays);
    SetLength(calcVj, n, NTrays);
    for j := 0 to NTrays-1 do
      begin
        calcTj[n-1, j] := Tj[j+1];
        calcLj[n-1, j] := Lj[j+1];
        calcVj[n-1, j] := Vj[j+1];
      end;
   
    for j := 1 to NTrays do
      begin
        Lj0[j] := calcLj[n-1, j-1];
        Vj0[j] := calcVj[n-1, j-1];
      end;
     
  until getErrorValue(n, calcTj, calcLj, calcVj) <= tolerance;
  ShowMessage(IntToStr(n))

  end;

end.
