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
  fu = function(temp, pres, Tc, Pc, om: double): double;

  TMatBalance = Class
    function TettaCorrection(bi__di: arrComp; F, D: double; xf: arrComp): double;
    procedure RelativeFugasity(T: double; var alpha: arrComp);
    procedure WilsonCorrelation(CritT, CritP, omega: arrComp; NTrays: integer;
      Tj, Pj: TArrOfDouble; var Kji, alpha: TArrOfArrOfDouble);
    function Wilson(Tj: double; Pj: double; Tc, Pc, om: double): double;
    function dihotomy(a, b: double; Pj: double; Tc, Pc, om: double): double;
    function getTj_0(Fj: arrTrays; zf: TArrOfArrOfDouble): arrTrays;
    function Kji_Recalc(Tj, Pj: TArrOfDouble; xji: TArrOfArrOfDouble): TArrOfDouble;
    function Tj_Recalc(Kj: TArrOfDouble; Tj, Pj: TArrOfDouble): TArrOfDouble;
    procedure OveralMatBalance(NTrays, FeedTray: integer; F, D, RefluxRate, Tcond,
      Treb, Pcond, Preb: double; xf: arrComp; var Res: TArrOfDouble);
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
    { Public declarations }
  End;

implementation

procedure TMatBalance.RelativeFugasity(T{K}: Double; var alpha: arrComp);
const
  a: arrComp = (31.35,	44.0103,	52.3785,	58.7845,	66.945,	66.7563,	63.3315,	70.4265,	78.3285);
  b: arrComp = (-1307.52,	-2568.82,	-3490.55,	-4136.68,	-4604.09,	-5059.18,	-5117.78,	-6055.6,	-6947);
  c: arrComp = (0,	0,	0,	0,	0,	0,	0,	0,	0);
  d: arrComp = (-3.26134,	-4.97635,	-6.10875,	-7.01666,	-8.25491,	-8.08935,	-7.48305,	-8.37865,	-9.44866);
  e: arrComp = (0.000029418,	0.0000146447,	0.0000111869,	0.0000103662,	0.0000115706,	0.00000925395,	0.00000776606,
                0.00000661666,	0.00000647481);
  f: arrComp = (2,	2,	2,	2,	2,	2,	2,	2,	2);
var
  i: integer;
  P: arrComp;
  min: double;

begin
  for i := 1 to NComp do
    P[i]:= a[i] + b[i] / (T + c[i]) + d[i] * ln(T) + e[i] * exp(f[i] * ln(T));
  min:= P[1];
  for i := 2 to NComp do
    if min > P[i] then
      min:= P[i];
  for i := 1 to NComp do
    alpha[i]:= P[i] / min;
end;

procedure TMatBalance.WilsonCorrelation(CritT: arrComp; CritP: arrComp; omega: arrComp;
  NTrays: Integer; Tj: TArrOfDouble; Pj: TArrOfDouble; var Kji, alpha: TArrOfArrOfDouble);
var
  i, j, k: integer;
  Ps: TArrOfArrOfDouble;
  min: double;

begin
  SetLength(Ps, NTrays+2, NComp);
  for j := 0 to NTrays+1 do
    begin
      min := 1e6;
      for i := 1 to NComp do
        begin
          Ps[j, i-1] := CritP[i] / 100 * exp(5.372697 * (1 + omega[i]) * (1 - (CritT[i] + 273.15) / Tj[j]));
          Kji[j, i-1] := Ps[j, i-1] {/ 100} / Pj[j];
          if min > Ps[j, i-1] then
            min := Ps[j, i-1];
        end;
      for k := 1 to NComp do
        alpha[j, k-1] := Ps[j, k-1] / min;
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

function TMatBalance.dihotomy(a: Double; b: Double; Pj: Double; Tc, Pc, om: double): double;
const
  eps = 1e-5;
var
  temp: double;
  function f(temp, pres, Tc, Pc, om: double): double;
  var
    i: integer;
    Ki: double;
  begin
    //s1 := 0;
    Ki := Wilson(temp, pres, Tc, Pc, om);
    Result := Ki - 1;
  end;
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
      Ti_sat[i] := dihotomy(90, 900, 0.195, Tcc[i], Pcc[i], omega[i]);
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

function TMatBalance.TettaCorrection(bi__di: arrComp; F, D: Double; xf: arrComp): double;
const
  eps = 1e-3;
var
  i: integer;
  x0: double;
  n: integer;
  function g(x: double): double;
  var
    i: integer;
  begin
    Result := 0;
    for i := 1 to NComp do
      Result := Result + F * xf[i] / (1 + x * (bi__di[i]));
    Result := Result - D;
  end;
  function g1(x: double): double;
  var
    i: integer;
  begin
    Result := 0;
    for i := 1 to NComp do
      Result := Result - bi__di[i] * F * xf[i] / sqr(1 + x * (bi__di[i]));
  end;
begin
  Result := 0;
  n := 0;
  repeat
    n:= n + 1;
    x0 := Result;
    Result := x0 - g(x0) / g1(x0);
    if n >= 20e3 then
      begin
        ShowMessage('Tetta-Method - Can not solve!');
        break
      end;
  until abs(Result - x0) <= eps;
end;

function TMatBalance.Kji_Recalc(Tj, Pj: TArrOfDouble; xji: TArrOfArrOfDouble): TArrOfDouble;
var
  i, j, k: integer;
  alpha: TArrOfArrOfDouble;
  s: double;
  Kj: TArrOfDouble;
  kji: TArrOfArrOfDouble;

begin
  SetLength(Kj, NTrays+2);
  SetLength(Kji, NTrays+2, NComp);
  SetLength(alpha, NTrays+2, NComp);
  for j := 0 to NTrays+1 do
    begin
      s:= 0;
      //RelativeFugasity(Tj[j], alpha);
      WilsonCorrelation(Tcc, Pcc, omega, NTrays, Tj, Pj, Kji, alpha);
      for i := 1 to NComp do
        s := s + xji[j, i-1] * alpha[j, i-1];
      Kj[j] := 1 / s;
      Result[j] := Kj[j]
      {
      for k := 1 to NComp do
        Result[j, k-1] := alpha[j, k-1] * Kj[j];   }
    end;
end;

function TMatBalance.Tj_Recalc(Kj: TArrOfDouble; Tj, Pj: TArrOfDouble): TArrOfDouble;
const
  T0 = 100;
  Tk = 900;
  h =  2.7315;
var
  i, j, k: integer;
  rTj: TArrOfDouble;
  deltaj: TArrOfDouble;
  deltaji: TArrOfArrOfDouble;
  rKji: TArrOfArrOfDouble;
  alpha: TArrOfArrOfDouble; // ������������� ���������
  temp: TArrOfDouble; // ���������� - �����
  err: double;
  s: double;
begin
  SetLength(rKji, NTrays+2, NComp);
  SetLength(deltaji, NTrays+2, NComp);
  SetLength(alpha, NTrays+2, NComp);
  SetLength(deltaj, NTrays+2);
  SetLength(temp, NTrays+2);
  SetLength(rTj, NTrays+2);
  for j := 0 to NTrays+1 do
    begin
      for i := 0 to NComp-1 do
        deltaji[j, i] := 0;
      rTj[j] := T0;
    end;
  k := 1;
  WilsonCorrelation(Tcc, Pcc, omega, NTrays, rTj, Pj, rKji, alpha);
  for j := 0 to NTrays+1 do
    begin
      deltaj[j] := 0;
      for i := 0 to NComp-1 do
        begin
          deltaji[j, i] := rKji[j, i] / alpha[j, i];
          deltaj[j] := {deltaj[j] +} deltaji[j, i] - Kj[j];
        end;
      temp[j] := deltaj[j];
    end;

  for j := 0 to NTrays+1 do
    begin
      deltaj[j] := 0;

      while (rTj[j] >= T0) and (rTj[j] <= Tk) do
        begin
          s := 0;
          rTj[j] := rTj[j] + h * k;
          WilsonCorrelation(Tcc, Pcc, omega, NTrays, rTj, Pj, rKji, alpha);
          for i := 0 to NComp-1 do
            begin
              deltaji[j, i] := rKji[j, i] / alpha[j, i];
              deltaj[j] := {deltaji[j, i]} deltaji[j, i] - Kj[j];
            end;
          s := s + (Kj[j] + deltaji[j, 1]) / Kj[j];
          if temp[j] > deltaj[j] then
            begin
              temp[j] := deltaj[j];
              Result[j] := rTj[j];
            end
          else
            Result[j] := Tj[j];
          err := sqrt(s) / (NTrays + 2);
        end;

    end;
end;

procedure TMatBalance.OveralMatBalance(NTrays, FeedTray: integer; F: Double; D: Double; RefluxRate: Double; Tcond: Double;
  Treb: Double; Pcond: Double; Preb: Double; xf: arrComp; var Res: TArrOfDouble);
const

  PhaseConst_Error = 1000;
var
  i: integer; // ���������
  j: integer; // �������
  k: integer; // iteration
  Tj, Pj: TArrOfDouble;
  Lj, Vj: TArrOfDouble;
  Kji: TArrOfArrOfDouble;
  Aji: TArrOfArrOfDouble;
  Sji: TArrOfArrOfDouble;
  xD: arrComp;
  xB: arrComp;
  di: arrComp;
  bi: arrComp;
  RefluxRatio: double;
  B: double;
  vji_di: TArrOfArrOfDouble;
  lji_di: TArrOfArrOfDouble;
  lji_bi: TArrOfArrOfDouble;
  vji_bi: TArrOfArrOfDouble;
  bi__di: arrComp;
  lji: TArrOfArrOfDouble;
  vji: TArrOfArrOfDouble;
  sd, sb, sl, sv: double;
  xji: TArrOfArrOfDouble;
  yji: TArrOfArrOfDouble;
  RecalcKji: TArrOfArrOfDouble;
  Error_Kji: double;
  rError_Kji: double;
  RecalcTj: TArrOfDouble;
  alpha: TArrOfArrOfDouble;
  Kj: TArrOfDouble;

  function getE(Kji, rKji: TArrOfArrOfDouble): double;
  var
    i, j: integer;
    s: double;
  begin
    s := 0;
    for j := 0 to NTrays+1 do
      for i := 0 to NComp-1 do
        s := s + (rKji[j, i] + Kji[j, i]) / Kji[j, i];
    Result := sqrt(s) / (NComp * (NTrays+2))
  end;

begin
  SetLength(Lj, NTrays+2);
  SetLength(Vj, NTrays+2);
  SetLength(Tj, NTrays+2);
  SetLength(Pj, NTrays+2);
  SetLength(Tj, NTrays+2);
  SetLength(Pj, NTrays+2);
  SetLength(Kji, NTrays+2, NComp);
  SetLength(Sji, NTrays+2, NComp);
  SetLength(Aji, NTrays+2, NComp);
  SetLength(vji_di, NTrays+2, NComp);
  SetLength(lji_di, NTrays+2, NComp);
  SetLength(vji_bi, NTrays+2, NComp);
  SetLength(lji_bi, NTrays+2, NComp);
  SetLength(vji, NTrays+2, NComp);
  SetLength(lji, NTrays+2, NComp);
  SetLength(xji, NTrays+2, NComp);
  SetLength(yji, NTrays+2, NComp);
  SetLength(RecalcKji, NTrays+2, NComp);
  SetLength(RecalcTj, NTrays+2);
  SetLength(alpha, NTrays+2, NComp);
  SetLength(Kj, NTrays+2);
  k := 0;
  Lj[0] := RefluxRate;
  Lj[NTrays+1] := F - D;
  for j := 2 to NTrays+1 do
    if j <> FeedTray then
      Lj[j-1] := Lj[j-2]
    else
      Lj[j-1] := Lj[j-2] + F;
  Vj[0] := 0; // ������ �����������
  Vj[1] := D + RefluxRate;
  for j := 2 to NTrays+1 do
    Vj[j]:= Vj[j-1];
  Tj[0] := Tcond;
  Tj[NTrays+1] := TReb;
  Pj[0] := Pcond;
  Pj[NTrays+1] := Preb;
  for j := 1 to NTrays do
    begin
      Tj[j] := Tj[j-1] + (Tj[NTrays+1] - Tj[0]) / (NTrays + 1);
      Pj[j] := Pj[j-1] + (Pj[NTrays+1] - Pj[0]) / (NTrays + 1);
    end;
  RefluxRatio := Lj[0] / D;
  B := F - D;

  for j := 0 to NTrays+1 do
    begin
      Tj[j] := Tj[j] + 273.15;
      Pj[j] := Pj[j] * 10;
    end;

  Repeat
  k := k + 1;
  WilsonCorrelation(Tcc, Pcc, omega, Ntrays, Tj, Pj, Kji, alpha);
 // ������ ����������� ������
  for j := 1 to FeedTray-2 do
    for i := 1 to NComp do
      Aji[j, i-1] := Lj[j] / (Kji[j, i-1] * Vj[j]);

  for i := 1 to NComp do
    vji_di[1, i-1] := Lj[0] / D + 1;

  for j := 2 to FeedTray-1 do
    for i := 1 to NComp do
      begin
        lji_di[j-1, i-1] := Aji[j-1, i] * vji_di[j-1, i-1];
        vji_di[j, i-1] := lji_di[j-1, i-1] + 1;
      end;
  // ������ ������������� ������
  for i := 1 to NComp do
    begin
      Sji[NTrays+1, i-1] := Kji[NTrays+1, i-1] * Vj[NTrays+1] / B;
      lji_bi[NTrays, i-1] := Sji[NTrays+1, i-1] + 1;
    end;
  for j := NTrays downto FeedTray-1 do
    for i := 1 to NComp do
      begin
        Sji[j, i-1] := Kji[j, i-1] * Vj[j] / Lj[j];
        vji_bi[j, i-1] := Sji[j, i-1] * lji_bi[j, i-1];
        lji_bi[j-1, i-1] := vji_bi[j, i-1] + 1;
      end;
  for i := 1 to NComp do
    begin
      bi__di[i] := vji_di[FeedTray-1, i-1] / vji_bi[FeedTray-1, i-1];
      di[i] := F * xf[i] / (1 + lji_di[0, i]);
    end;
  //ShowMessage('tetta = ' + FloatToStr(RoundTo(TettaCorrection(bi__di, F, D, xf), -8)));
  sd := 0;
  sb := 0;
  for i := 1 to NComp do
    begin
      di[i] := F * xf[i] / (1 + TettaCorrection(bi__di, F, D, xf) * bi__di[i]);
      sd := sd + di[i];
      bi[i] := TettaCorrection(bi__di, F, D, xf) * bi__di[i] * di[i];
      sb := sb + bi[i];
    end;
  D := sd;
  B := sb;
  // ������ ����� �������� ������� �� ���� � ��������
  for j := 1 to FeedTray-1 do
    begin
      sl := 0;
      sv := 0;
      for i := 1 to NComp do
        begin
          vji[j, i-1] := vji_di[j, i-1] * di[i];
          sv := sv + vji[j, i-1];
          lji[j, i-1] := lji_di[j, i-1] * di[i];
          sl := sl + lji[j, i-1];
        end;
      Vj[j] := sv;
      Lj[j] := sl;
    end;
  for j := FeedTray-1 to NTrays do
    begin
      sl := 0;
      sv := 0;
      for i := 1 to NComp do
        begin
          vji[j, i-1] := vji_bi[j, i-1] * bi[i];
          sv := sv + vji[j, i-1];
          lji[j, i-1] := lji_bi[j, i-1] * bi[i];
          sl := sl + lji[j, i-1];
        end;
      Vj[j] := sv;
      Lj[j] := sl;
    end;
  // ������ ��������
  for i := 1 to NComp do
    begin
      xD[i] := di[i] / D;
      xB[i] := bi[i] / B;
      for j := 1 to NTrays+1 do
        begin
          xji[j, i-1] := lji[j, i-1] / Lj[j];
          yji[j, i-1] := vji[j, i-1] / Vj[j];
        end;
      xji[0, i-1] := xD[i];
      xji[NTrays+1, i-1] := xB[i];
    end;

   // ������ ����� �������� Kji
  {RecalcKji}
  Kj := Kji_Recalc(Tj, Pj, xji);
  RecalcTj := Tj_Recalc({RecalcKji}Kj, Tj, Pj);
  for j := 0 to NTrays+1 do
    Tj[j] := RecalcTj[j];

  for j := 0 to NTrays+1 do
    for i := 0 to NComp-1 do
      RecalcKji[j, i] := alpha[j, i] * Kj[j];

  rError_Kji := Error_Kji;
  Error_Kji := 0;
  for j := 0 to NTrays+1 do
    for i := 0 to NComp-1 do
      Error_Kji := Error_Kji + abs(Kji[j, i] - RecalcKji[j, i]);

  //ShowMessage(FloatToStr(getE(Kji, RecalcKji)));

  for j := 0 to NTrays+1 do
    for i := 0 to NComp-1 do
      if k <= 5 then
        Kji[j, i] := RecalcKji[j, i]
      else
        begin
          Kji[j, i] := sqrt(Kji[j, i] * RecalcKji[j, i]);
          k := 1;
        end;

  until {Error_Kji <= PhaseConst_Error} abs(rError_Kji - Error_Kji) <= 1 {getE(Kji, RecalcKji) <= 1};

end;

procedure TMatBalance.MatBalCalculation(Fl: arrTrays; Fv: arrTrays; Wl: arrTrays;
  Wv: arrTrays; T1: Double; TN: Double; P1: Double; PN: Double; var L: arrTrays; var V: arrTrays);
var
  zf: TArrOfArrOfDouble;
  Res: TArrOfDouble;
  Fj: arrTrays;
  Tj_0: arrTrays;
  i, j: Integer;
begin
  SetLength(zf, NComp, NTrays);

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
    Fj[j] := 0;
  Fj[FeedTray] := 13.8;

  //OveralMatBalance(10, 5, 13.8, 4.5, 13.5, 60, 100, 0.19, 0.20, zf, Res);

  Tj_0 := getTj_0(Fj, zf);
end;

end.
