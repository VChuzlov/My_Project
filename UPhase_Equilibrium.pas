unit UPhase_Equilibrium;

interface

const
  NComp = 9;


type
  arrComp = array[1..NComp] of double;


  TPhase_Eq_Calc = class

    //procedure WilsonCorrelation(CritT, CritP, omega: arrComp; T, P: double; var K: arrComp);
    procedure PhaseEqCalculation;

    private
      { Private declarations }
    public
      CritT: arrComp;
      CritP: arrComp;
      omega: arrComp;
      T: Double;
      P: Double;
      K: arrComp;
      { Public declarations }
  end;

implementation

 {
procedure TPhase_Eq_Calc.WilsonCorrelation(CritT: arrComp; CritP: arrComp;
  omega: arrComp; T: Double; P: Double; var K: arrComp);
const
  h = 1e-5;
  eps = 1e-3;
var
  i: integer;
  Ps: arrComp; // Парциальное давление пара
  sl, sv: double; // Сумма мольных долей жидкости и пара
  n: integer; // число итераций
  a, b: double;
begin
  for i := 1 to NComp do
      begin
        CritP[i]:= random;
        CritT[i]:= random;
        omega[i]:= random;
      end;
  P := P * 10; // Перевод в бары
  T := T + 273.15;
  for i := 1 to NComp do
    begin
      Ps[i]:= CritP[i] * exp(5.372697 * (1 + omega[i]) * (1 - CritT[i] / T));
      K[i]:= Ps[i] / P;
    end;
end;
         }
procedure TPhase_Eq_Calc.PhaseEqCalculation;
begin
  //
end;

end.
