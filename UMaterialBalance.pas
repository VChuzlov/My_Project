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
    procedure WilsonCorrelation(CritT, CritP, omega: arrComp; NTrays: integer;
      Tj, Pj: arrTrays; var Kij: TArrOfArrOfDouble);
    function Wilson(Tj: double; Pj: double; Tc, Pc, om: double): double;
    function dihotomy(a, b: double; Pj: double; Tc, Pc, om: double): double;
    function getTj_0(Fj: arrTrays; zf: TArrOfArrOfDouble): arrTrays;
    procedure InitGuessLV(Fj, Uj, Wj: arrTrays; WD, LD, L0: double; var Lj0, Vj0: arrTrays);
    procedure Gauss_Jordan(arr: TArrOfArrOfDouble; var x: TArrOfDouble);
    procedure CalculateLiquidMoleFractions(Fj, Lj, Vj, Uj, Wj: arrTrays; zij, Kij: TArrOfArrOfDouble;
      var xij: TArrOfArrOfDouble);
    procedure NormalizeOfxij(xij: TArrOfArrOfDouble; var norm_xij: TArrOfArrOfDouble);
    procedure Secant(Tj: arrTrays; xij, Kij: TArrOfArrOfDouble; var rTj: arrTrays);
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

procedure TMatBalance.InitGuessLV(Fj: arrTrays; Uj: arrTrays; Wj: arrTrays;
  WD: Double; LD: Double; L0: Double; var Lj0: arrTrays; var Vj0: arrTrays);
var
  i, j: integer;
  rD, rB: double;
  dj: arrTrays;
  qj: arrTrays; // feed quality
begin
  Lj0[1] := Fj[1] + L0 - Uj[1];
  Vj0[1] := 0;
  for j := 2 to NTrays do
    Lj0[j] := Lj0[j-1] + Fj[j] - Uj[j];
  for j := 1 to NTrays do
    Vj0[j] := WD + LD + L0 - Wj[j];
  for j := 1 to NTrays do
    begin
      dj[j] := (Uj[j] + Wj[j]) / (Vj0[j] + Lj0[j]);
      qj[j] := 1;
    end;
  for j := NTrays-1 downto 2 do
    Vj0[j] := ((1 - qj[j]) * Fj[j] + Vj0[j+1]) / (dj[j] + 1);
  for j := 2 to NTrays-1 do
    Lj0[j] := (Fj[j] + Vj0[j+1] + Lj0[j-1]) / (dj[j] + 1) - Vj0[j];
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

procedure TMatBalance.NormalizeOfxij(xij: TArrOfArrOfDouble; var norm_xij: TArrOfArrOfDouble);
var
  i, j: integer;
  s: TArrOfDouble;
begin
  SetLength(s, NTrays);
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
        norm_xij[i, j] := xij[i, j] / s[j];
  end;

procedure TMatBalance.Secant(Tj: arrTrays; xij, Kij: TArrOfArrOfDouble; var rTj: arrTrays);
const
  eps = 1e-5;
  function f(T: double): double;
    var
      i: integer;
      s: double;
  begin
    s := 0;

  end;
begin
  //
end;

procedure TMatBalance.MatBalCalculation(Fl: arrTrays; Fv: arrTrays; Wl: arrTrays;
  Wv: arrTrays; T1: Double; TN: Double; P1: Double; PN: Double; var L: arrTrays; var V: arrTrays);
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

begin
  SetLength(zf, NComp, NTrays);
  SetLength(Kij, NComp, NTrays);
  SetLength(xij, NComp, NTrays);
  SetLength(norm_xij, NComp, NTrays);
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
    end;
  Fj[FeedTray] := 13.8;
  Pj[1] := P1;
  Pj[Ntrays] := PN;
  for j := 2 to NTrays-1 do
    Pj[j] := Pj[j-1] + (Pj[NTrays] - Pj[1]) / NTrays;

  Tj_0 := getTj_0(Fj, zf);
  InitGuessLV(Fj, Uj, Wj, 0, 4.5, 13.5, Lj0, Vj0);
  WilsonCorrelation(Tcc, Pcc, omega, NTrays, Tj_0, Pj, Kij);
  CalculateLiquidMoleFractions(Fj, Lj0, Vj0, Uj, Wj, zf, Kij, xij);
  NormalizeOfxij(xij, norm_xij);

end;

end.
