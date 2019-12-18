unit UMaterialBalance;

interface
uses
  UPhase_Equilibrium;

const
  NTrays = 12;
  FeedTray = 5;

type
  arrTrays = array[1..NTrays] of double;

  TMatBalance = Class
    procedure RelativeFugasity(T: double; var alpha: arrComp);
    procedure OveralMatBalance(Fl, Fv, Wl, Wv: arrTrays; var L, V: arrTrays);
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

procedure TMatBalance.OveralMatBalance(Fl: arrTrays; Fv: arrTrays; Wl: arrTrays;
  Wv: arrTrays; var L: arrTrays; var V: arrTrays);
var
  i: integer;
  dV: arrTrays; // V[i] - V[i+1]
begin
  L[1] := 1000;
  L[NTrays] := 900;
  for i := 2 to NTrays-1 do
    L[i] := L[i-1] + (L[1] - L[NTrays]) / NTrays;
  for i := 2 to NTrays-1 do
    dV[i] := L[i-1] - L[i] + Fl[i] + Fv[i] - Wl[i-1] - Wv[i+1];
  dV[1] :=  Fl[1] + Fv[1] - Wv[2] - L[1];
  dV[NTrays] := L[NTrays-1] - Wl[NTrays-1] + Fl[NTrays] + Fv[NTrays] - L[NTrays];
  V[NTrays] := dV[NTrays];
  for i := NTrays-1 downto 1 do
    V[i] := dV[i] + dV[i+1];
end;

procedure TMatBalance.MatBalCalculation(Fl: arrTrays; Fv: arrTrays; Wl: arrTrays;
  Wv: arrTrays; T1: Double; TN: Double; P1: Double; PN: Double; var L: arrTrays; var V: arrTrays);
var
  i: integer;
  Tprofile: arrTrays;
  Pprofile: arrTrays;
  Kij: array [1..NTrays, 1..NComp] of double;
  j: Integer;
  PhaseCalc: TPhase_Eq_Calc;
begin
  for i := 1 to NTrays do
    begin
      Fl[i] := 0;
      Fv[i] := 0;
      Wl[i] := 0;
      Wv[i] := 0;
    end;
  Fl[FeedTray] := 1000;
  Wl[1]:= 500;
  Wl[NTrays] := 500;
  OveralMatBalance(Fl, Fv, Wl, Wv, L, V);
  TProfile[1] := T1;
  TProfile[NTrays] := TN;
  PProfile[1] := P1;
  PProfile[NTrays] := PN;
  for i := 2 to NTrays-1 do
    begin
      TProfile[i]:= TProfile[i-1] + (TProfile[1] - TProfile[NTrays]) / NTrays;
      PProfile[i]:= PProfile[i-1] + (PProfile[1] - PProfile[NTrays]) / NTrays;
    end;
  PhaseCalc:= TPhase_Eq_Calc.Create;

  for i := 1 to NTrays do
    for j := 1 to NComp do
      PhaseCalc.WilsonCorrelation(PhaseCalc.CritT, PhaseCalc.CritP, PhaseCalc.omega, TProfile[j], PProfile[j], PhaseCalc.K);

end;

end.
