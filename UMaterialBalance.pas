unit UMaterialBalance;

interface
uses
  UPhase_Equilibrium, Dialogs, System.SysUtils, Math;

const
  NTrays = 10;
  FeedTray = 5;

type
  arrTrays = array[1..NTrays] of double;
  TArrOfDouble = array of double;
  TArrOfArrOfDouble = array of array of double;

  TMatBalance = Class
    function TettaCorrection(bi__di: arrComp; F, D: double; xf: arrComp): double;
    procedure RelativeFugasity(T: double; var alpha: arrComp);
    procedure WilsonCorrelation(CritT, CritP, omega: arrComp; NTrays: integer;
      Tj, Pj: TArrOfDouble; var Kji: TArrOfArrOfDouble);
    function Kji_Recalc(Tj: TArrOfDouble; xji: TArrOfArrOfDouble): TArrOfArrOfDouble;
    procedure OveralMatBalance(NTrays, FeedTray: integer; F, D, RefluxRate, Tcond,
      Treb, Pcond, Preb: double; xf: arrComp; var Res: TArrOfDouble);
    procedure MatBalCalculation(Fl, Fv, Wl, Wv: arrTrays; T1, TN, P1, PN: double;
      var L, V: arrTrays);
  private
    { Private declarations }
  public
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
  NTrays: Integer; Tj: TArrOfDouble; Pj: TArrOfDouble; var Kji: TArrOfArrOfDouble);
var
  i, j: integer;
  Ps: TArrOfArrOfDouble;
begin
  SetLength(Ps, NTrays+2, NComp);

  for j := 0 to NTrays+1 do
    begin
      Tj[j] := Tj[j] + 273.15;
      Pj[j] := Pj[j] * 10;
    end;
  for j := 0 to NTrays+1 do
    for i := 1 to NComp do
      begin
        Ps[j, i-1] := CritP[i] * exp(5.372697 * (1 + omega[i]) * (1 - (CritT[i] + 273.15) / Tj[j]));
        Kji[j, i-1] := Ps[j, i-1] / 100 / Pj[j];
      end;
end;

function TMatBalance.TettaCorrection(bi__di: arrComp; F, D: Double; xf: arrComp): double;
const
  eps = 1e-5;
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

function TMatBalance.Kji_Recalc(Tj: TArrOfDouble; xji: TArrOfArrOfDouble): TArrOfArrOfDouble;
var
  i, j, k: integer;
  alpha: arrComp;
  s: double;
  Kj: TArrOfDouble;

begin
  SetLength(Kj, NTrays+2);
  for j := 0 to NTrays+1 do
    begin
      s:= 0;
      RelativeFugasity(Tj[j], alpha);
      for i := 1 to NComp do
        s := s + xji[j, i-1] * alpha[i];
      Kj[j] := 1 / s;
      for k := 1 to NComp do
        Result[j, k-1] := alpha[k] * Kj[j];
    end;
end;

procedure TMatBalance.OveralMatBalance(NTrays, FeedTray: integer; F: Double; D: Double; RefluxRate: Double; Tcond: Double;
  Treb: Double; Pcond: Double; Preb: Double; xf: arrComp; var Res: TArrOfDouble);
const
  Tcc: arrComp = (-82.45, 32.28, 96.75, 134.9, 152, 196.5, 187.2, 234.7, 267);
// ����������� ��������
  Pcc: arrComp = (4641, 4484, 4257, 3648, 3797, 3375, 3334, 3032, 2737);
// ������������� ������
  omega: arrComp = (0.0115, 0.0986, 0.1524, 0.1848, 0.201, 0.2539, 0.2222, 0.3007, 0.3498);
var
  i: integer; // ���������
  j: integer; // �������
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
  WilsonCorrelation(Tcc, Pcc, omega, Ntrays, Tj, Pj, Kji);
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
  RecalcKji := Kji_Recalc(Tj, xji);
end;

procedure TMatBalance.MatBalCalculation(Fl: arrTrays; Fv: arrTrays; Wl: arrTrays;
  Wv: arrTrays; T1: Double; TN: Double; P1: Double; PN: Double; var L: arrTrays; var V: arrTrays);
var
  xf: arrComp;
  Res: TArrOfDouble;

begin
  xf[1] := 2.22222222222222e-002;
  xf[2] := 4.44444444444444e-002;
  xf[3] := 6.66666666666667e-002;
  xf[4] := 8.88888888888889e-002;
  xf[5] := 0.111111111111111;
  xf[6] := 0.133333333333333;
  xf[7] := 0.155555555555556;
  xf[8] := 0.177777777777778;
  xf[9] := 0.200000000000000;
  OveralMatBalance(10, 5, 13.8, 4.5, 13.5, 60, 100, 0.19, 0.20, xf, Res);
end;

end.
