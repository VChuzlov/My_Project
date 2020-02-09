unit UMaterialBalance;

interface
uses
  UPhase_Equilibrium, Dialogs, System.SysUtils, Math;

const
  NTrays = 62;
  FeedTray1 = 10;
  FeedTray2 = 22;

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
      var Uj: arrTrays; var Wj: arrTrays; var dj: arrTrays; var rB: double);
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
      Qj: arrTrays; ej: arrTrays; LD, VD: double; rB: double; var Vj: arrTrays);
    procedure CalculateHeatDuties (Fj, ej, Lj, Vj, Uj, Wj, Hf_l, Hf_v, H_l, H_v: arrTrays;
      var Qc, Qr: double);
    procedure CalculateSideFlows_and_LiquidFlows(dj, Fj, ej, Vj, Lj0: arrTrays; LD, L0: double; zf, xij: TArrOfArrOfDouble;
      var Wj, Uj, Lj: arrTrays);
    function getErrorValue(n: integer; calcTj, calcLj, CalcVj: TArrOfArrOfDouble): double;

    function forRegularTrays(T, P: double; xi: TArrOfDouble): double;
    function forCondenser (T, P: double; zf: TArrOfDouble): double;
    function teta_method(Fj: arrTrays; zf, xij: TArrOfArrOfDouble; D, B, Dco: double): double;
    procedure LV_correction(rD, rB: double; Fj, Uj, Wj: arrTrays; var Lj, Vj: arrTrays);
    procedure correction_with_tray_efficiencies(tray_efficiency: arrTrays; Kij: TarrOfArrOfDouble; var xij, yij: TArrOfArrOfDouble);
    procedure MatBalCalculation(Fl, Fv, Wl, Wv: arrTrays; T1, TN, P1, PN: double;
      D: double; LD: double;
      var Tj: arrTrays; var Lj: arrTrays; var Vj: arrTrays;
      var xij: TArrOfArrOfDouble; var yij: TArrOfArrOfDouble;
      var calcTj: TArrOfArrOfDouble;
      var calcLj: TArrOfArrOfDouble;
      var calcVj: TArrOfArrOfDouble;
      var n: integer);
  private
    { Private declarations }
  public
    const
      // ����������� �����������
      Tcc: arrComp = ({-82.45, 32.28, 96.75, 134.9, 152, 187.2, 196.5, 234.7, 267}
                      146.450006103516,	155.477014160156,	144.748010253906,
                      134.946008300781,	152.049005126953,	270.810021972656,
                      300.409020996094,	293.258020019531,	276.909020996094,
                      280.370019531250,	290.339013671875,	91.8500000000000,
                      246.639001464844,	258.019006347656,	257.219018554688,
                      264.198022460938,	262.100000000000,	247.350000000000,
                      286.488000488281,	191.549005126953,	298.198022460938,
                      91.8500000000000,	187.248010253906,	226.830010986328,
                      215.850000000000,	224.347009277344,	231.299005126953,
                      385.149011230469,	337.222009277344,	365.149011230469,
                      32.2780090332031
                      );
      // ����������� ��������
      Pcc: arrComp = ({4641, 4484, 4257, 3648, 3797, 3375, 3334, 3032, 2737}
                      4022.60009765625,	4102.35009765625,	4002.33007812500,
                      3647.62011718750,	3796.62011718750,	2567.57006835938,
                      2819.87011718750,	2729.62011718750,	2486.51000976563,
                      2556.37011718750,	2628.37011718750,	4620.41015625000,
                      2736.78002929688,	2953.62011718750,	2733.62011718750,
                      2908.02001953125,	2813.79003906250,	2773.26000976563,
                      2484.37011718750,	3528.69995117188,	2479.37011718750,
                      4256.66015625000,	3333.59008789063,	3126.87011718750,
                      3100.00000000000,	3010.36010742188,	3123.84008789063,
                      1829.92004394531,	2095.87011718750,	1964.93005371094,
                      4883.85009765625
                      );
      // ������������� ������
      omega: arrComp = ({0.0115, 0.0986, 0.1524, 0.1848, 0.201, 0.2539, 0.2222, 0.3007, 0.3498}
                        0.187000006437302,	0.210565000772476,	0.189980000257492,
                        0.184790000319481,	0.20100000500679,	0.310000002384186,
                        0.28999000787735,	0.319990009069443,	0.345990002155304,
                        0.340990006923676,	0.340000003576279,	0.148000001907349,
                        0.307000011205673,	0.259990006685257,	0.340000003576279,
                        0.305000007152557,	0.326990008354187,	0.314000010490417,
                        0.38400000333786,	0.232960000634193,	0.416680008172989,
                        0.152400001883507,	0.222240000963211,	0.246950000524521,
                        0.231940001249313,	0.279100000858307,	0.275000005960464,
                        0.561990022659302,	0.409990012645721,	0.535000026226044,
                        0.0986000001430511
                        );

     // Ideal Gas Comp Heat Capasity
     IdealGasCp_a: arrComp = ({4.598,	1.292,	-1.009,	-0.332,	2.266,	-2.275,	-0.866,	-1.054,	-1.229}
                              -0.715,	0.105,	3.834,	-0.332,	2.266,	-2.201,
                              -2.201,	-2.201,	-2.201,	-2.201,	-2.201,	0.886,
                              -11.966,	-5.48,	-9.408,	-11.966,	-1.683,
                              -9.408,	-2.201,	-0.032,	-2.201,	-1.009,	-2.275,
                              -3.489,	-3.973,	-2.524,	0.57,	-14.932,	-3.928,
                              -7.473,	1.292
                              );
     IdealGasCp_b: arrComp = ({0.01245,	0.04254,	0.07315,	0.09189,	0.07913,	0.121,	0.1164,	0.139,	0.1615}
                              0.08436,	0.07054,	0.06698,	0.09189,	0.07913,
                              0.1877,	0.1877,	0.1877,	0.1877,	0.1877,	0.1877,
                              0.05602,	0.2139,	0.1796,	0.2064,	0.2139,	0.1633,
                              0.2064,	0.1877,	0.1034,	0.1877,	0.07315,	0.121,
                              0.1469,	0.1503,	0.1477,	0.1359,	0.2362,	0.1671,
                              0.1788,	0.04254
                              );
     IdealGasCp_c: arrComp = ({0.00000286,	-0.00001657,	-0.00003789,	-0.00004409,	-0.00002647,	-0.00006519,
                              -0.00006163,	-0.00007449,	-0.0000872}
                              -0.00004754,	-0.00002431,	-0.00002607,	-0.00004409,
                              -0.00002647,	-0.0001051,	-0.0001051,	-0.0001051,
                              -0.0001051,	-0.0001051,	-0.0001051,	-0.00002771,
                              -0.0001519,	-0.0001056,	-0.0001502,	-0.0001519,
                              -0.00008919,	-0.0001502,	-0.0001051,	-0.00005534,
                              -0.0001051,	-0.00003789,	-0.00006519,	-0.00008063,
                              -0.00008314,	-0.00008533,	-0.00006854,	-0.0001384,
                              -0.00009841,	-0.0001099,	-0.00001657
                              );
     IdealGasCp_d: arrComp = ({-0.000000002703,	0.000000002081,	0.000000007678,	0.000000006915,	-0.000000000647,
                              0.00000001367,	0.00000001267,	0.00000001551,	0.00000001829}
                              0.00000001066,	-0.000000000147,	0.000000002173,
                              0.000000006915,	-0.000000000674,	0.00000002316,
                              0.00000002316,	0.00000002316,	0.00000002316,
                              0.00000002316,	0.00000002316,	0.000000005266,
                              0.00000004146,	0.000000024,	0.00000004386,
                              0.00000004146,	0.00000001871,	0.00000004386,
                              0.00000002316,	0.00000001118,	0.00000002316,
                              0.000000007678,	0.00000001367,	0.00000001629,
                              0.00000001636,	0.00000001931,	0.00000001202,
                              0.00000003084,	0.00000002228,	0.00000002582,
                              0.000000002081
                              );

     // Liquid Comp Heat Capasity
     LiquidCp_a: array [1..8] of double = (10.1273,	-15.3546,	3.2008,	19.7302,	-0.8949,	-0.01489,	0.2241,	-0.04342);
     LiquidCp_b: array [1..5] of double = (0.31446,	2.5346,	-2.0242,	-0.07055,	0.07264);
     LiquidCp_R: arrComp = ({1.124,	1.8314,	2.4255,	2.8962,	2.8885,	3.313,	3.385,	3.812,	4.2665}
                            2.7458,	2.7765,	2.8281,	2.8962,	2.8885,	4.1714,
                            4.0859,	4.2052,	4.5932,	4.3463,	4.4084,	2.2283,
                            3.9634,	3.696,	4.2779,	3.921,	4.1454,	3.7952,
                            4.7401,	3.1956,	4.1556,	2.4255,	3.313,	3.5209,
                            3.4846,	3.809,	3.6797,	6.4321,	5.539,	5.9867,
                            1.8314
                            );
     LiquidCp_k: arrComp = ({0,	0,	0,	-0.6884,	0,	-0.7643,	0,	0,	0}
                            -0.0103,	0.2014,	-0.1428,	-0.6884,	0,	-1.3313,
                            -1.4284,	-1.0783,	-1.1169,	-0.7135,	-0.8089,
                            0.359,	-0.6858,	-1.3293,	-0.8844,	-0.7324,
                            -0.6597,	-1.161,	-0.9121,	1.1055,	-1.893,	0,
                            -0.7643,	-0.8695,	-1.1611,	-0.893,	-0.6343,
                            0,	0,	0,	0
                            );

     // Heat of Formation [kJ / kmole]
     dHf298: arrComp = ({-74900,	-84738,	-103890,	-134590,	-126190,	-154590,	-146490,	-167290,	-187890}
                        -125.970001220703,	-11178.900390625,	-16909,	-134590,
                        -126190,	-224290,	-216590,	-217590,	-222790,	-219590,
                        -214089,	20429,	-202090,	-204890,	-195090,	-199390,
                        -192390,	-206290,	-215590,	-20929,	-242090,	-103890,
                        -154590,	-177890,	-185690,	-174390,	-171690,	-291090,
                        -256510,	-270490,	-84738
                        );

     // Molar weight
     Mr: arrComp = ({16.0429000854492,	30.0699005126953,	44.0970001220703,	58.1240005493164,
                    58.1240005493164,	72.1510009765625,	72.1510009765625,	86.1779022216797,	100.205001831055}
                    56.1077003479004,	56.1077003479004,	56.1077003479004,
                    58.1240005493164,	58.1240005493164,	114.232002258301,
                    114.232002258301,	114.232002258301,	114.232002258301,
                    114.232002258301,	114.232002258301,	42.0806007385254,
                    100.205001831055,	100.205001831055,	100.205001831055,
                    100.205001831055,	100.205001831055,	100.205001831055,
                    114.232002258301,	70.1350009765625,	128.259002685547,
                    44.0970001220703,	72.1510009765625,	86.1779022216797,
                    86.1779022216797,	86.1779022216797,	86.1779022216797,
                    170.339004516602,	142.285003662109,	156.313003540039,
                    30.0699005126953
                    );

    // Normal Boiling Point
    tbpi: arrComp = ({-161.525,	-88.5999969482422,	-42.1019958496094,	-11.7299865722656,
                     -0.501989746093727,	27.8780151367188,	36.0590148925781,	68.7300048828125,
                     98.4290100097656}
                     -6.252,	0.877008056640648,	-6.85098876953123,	-11.7299865722656,
                     -0.501989746093727,	99.2380004882813,	114.764001464844,	113.472009277344,
                     109.106011962891,	109.432000732422,	115.610009765625,	-47.7509979248047,
                     80.4930053710938,	80.8760009765625,	90.0490051269531,	89.7780090332031,
                     91.8470092773438,	79.1910034179688,	117.653009033203,	29.9460083007813,
                     122.291009521484,	-42.1019958496094,	27.8780151367188,	57.9770141601563,
                     49.7310119628906,	60.2610107421875,	63.2700134277344,	216.278009033203,
                     167.028009033203,	195.890008544922,	-88.5999969482422
                     );

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
      Ti_sat[i] := dihotomy(getTsat, 90, 900, 0.645, Tcc[i], Pcc[i], omega[i]);
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
  var Uj: arrTrays; var Wj: arrTrays; var dj: arrTrays; var rB: double);
var
  i, j: integer;
  {dj: arrTrays;}
  qj: arrTrays; // feed quality
  rD: double; // ��������� ����� // ������ �� 1
  {rB: double; // ������� ����� // ������ �� 1 }
  D, B: double;
  s: double;
  TrayMaterialBalanceError: arrTrays;
  teta: double;

  procedure Calculation(rB: double; var Fj, Uj, Wj, Lj0, Vj0, dj: arrTrays);
  var
    j: integer;

  begin
    //Lj0[1] := L0;
    Lj0[1] := rD * LD;
    for j := 2 to NTrays-1 do
    Lj0[j] := Lj0[j-1] + Fj[j] - Uj[j];
    Lj0[NTrays] := Uj[NTrays]; // ����� ��������

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
      Calculation(rB, Fj, Uj, Wj, Lj0, Vj0, dj);
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
  rB := {1.93548387096774}1;
  Uj[1] := LD;
  Wj[1] := WD;
  Uj[NTrays] := Fj[FeedTray1] + Fj[FeedTray2] - (LD + WD);  // �.�. ����� Fj
  Vj0[1] := 0; // ������ �������� �� ���� ������������
  Lj0[1] := LD * rD;

  Calculation(rB, Fj, Uj, Wj, Lj0, Vj0, dj);
  rB := get_rB(1e-5, 1000, Fj, Uj, Wj, Lj0, Vj0);
  Calculation(rB, Fj, Uj, Wj, Lj0, Vj0, dj);

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
        s := s + zf[i] / (Ki {- 1});
      end;
    Result := 1 - s {- 1};
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

    function get_borders(x: double): double;
    const
      h = 1;
    begin
      while f(x, P, xi) * f(x + h, P, xi) / abs(f(x + h, P, xi)) >= 0 do
        x := x + h;
      Result := x + h;
    end;

  begin
    b := get_borders(a);
    if f(a, P, xi) * f(b, P, xi) / abs(f(b, P, xi)) < 0 then
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
          DihotomyIterations(forCondenser, a, b, tmp);
        end
      else
        //SecantIterations(forRegularTrays, temp);
        DihotomyIterations(forRegularTrays, a, b, tmp);

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
  a: arrComp = ({-12.98,	-1.7675,	39.4889,	30.903,	67.721,	64.25,	63.198,
                74.513,	71.41}
                0.00000000830081,	0.0000000268367,	0.000000020158,	30.903,
                67.721,	122.94,	83.412,	77.083,	59.136,	63.746,	67.87,
                0.0000000192663,	75.9248,	86.0169,	47.7379,	80.006,
                0.0000000195712,	77.686,	125.477,	-0.0000000144437,
                67.8938,	39.4889,	64.25,	0,	0,	111.47,	83.821,
                0.0000000066763,	33.475,	0.00000000410848,	-1.7675
                );
  b: arrComp = ({2.36459,	1.1429,	0.395,	0.1533,	0.00854058,	-0.131798,
                -0.0117017,	-0.096697,	-0.0968949}
                -0.0533615,	0.3265,	0.28605,	0.1533,	0.00854058,	-0.364398,
                0.094612,	0.109999,	0.209496,	0.1845,	0.15289,	0.0881636,
                0.222689,	0.154798,	-0.125,	0.157998,	-0.056075,	0.215498,
                -0.060366,	-0.0019105,	0.14099,	0.395,	-0.131798,	-0.1695,
                -0.193,	-0.6057,	-0.1695,	-0.0547613,	0.209496,	-0.0537062,	1.1429
                );
  c: arrComp = ({-0.00213247,	-0.0003236,	0.00211409,	0.00263479,	0.00327699,
                0.003541,	0.0033164,	0.00347649,	0.003473}
                0.0031475,	0.00228487,	0.00249875,	0.00263479,	0.00327699,
                0.0042662,	0.0028405,	0.0028378,	0.0028279,	0.0028305,
                0.002835,	0.0027863,	0.00281928,	0.0028265,	0.00360259,
                0.0028285,	0.0033777,	0.0028195,	0.0034085,	0.00308619,
                0.0028379,	0.00211409,	0.003541,	0.00356839,	0.003651,
                0.0049203,	0.003678,	0.00337268,	0.0028405,	0.00337144,	-0.0003236
                );
  d: arrComp = ({0.0000056618,	0.0000042431,	0.000000396486,	0.0000000727226,
                -0.00000110968,	-0.0000013332,	-0.0000011705,	-0.0000013212,	-0.0000013302}
                -0.00000118222,	-0.000000416634,	-0.000000648153,	0.0000000727226,
                -0.00000110968,	-0.00000228148,	-0.000000686766,	-0.000000684396,
                -0.000000675566,	-0.000000677926,	-0.000000682036,	-0.000000918863,
                -0.000000667846,	-0.000000674316,	-0.0000012797,	-0.000000676066,
                -0.00000121066,	-0.000000668226,	-0.000001235,	-0.0000011012,
                -0.000000684406,	0.000000396486,	-0.0000013332,	-0.00000235029,
                -0.0000024234,	-0.000003017,	-0.0000015579,	-0.00000124201,
                -0.000000686766,	-0.00000123662,	0.0000042431
                );
  e: arrComp = ({-0.00000000372476,	-0.00000000339316,	-0.000000000667176,
                -0.000000000727896,	0.000000000176646,	0.000000000251446,
                0.000000000199636,	0.000000000252365,	0.000000000255766}
                0.000000000198854,	-0.0000000000400519,	0.0000000000405376,
                -0.000000000727896,	0.000000000176646,	0.000000000859866,
                0,	0,	0,	0,	0,	0.000000000130992,	0,	0,	0.0000000000986336,
                0,	0.000000000184802,	0,	0.000000000172896,	0.000000000166852,
                0,	-0.000000000667176,	0.000000000251446,	0.000000000641015,
                0.000000000643776,	0.00000000106506,	0.000000000353795,
                0.000000000199451,	0,	0.000000000197836,	-0.00000000339316
                );
  f: arrComp = ({0.000000000000860896,	0.000000000000882096,	0.000000000000167936,
                0.000000000000236736,	-6.39926E-15,	-1.29576E-14,	-8.66485E-15,
                -1.34666E-14,	-1.37726E-14}
                -4.27421E-23,	-1.48189E-22,	-1.09334E-22,	0.000000000000236736,
                -6.39926E-15,	-0.000000000000156206,	0,	0,	0,	0,	0,
                -1.04978E-22,	0,	0,	2.31836E-14,	0,	-1.08753E-22,	0,
                8.06576E-15,	8.01642E-23,	0,	0.000000000000167936,
                -1.29576E-14,	0,	0,	-0.000000000000108026,	-1.77176E-14,
                -3.4102E-23,	0,	-2.08778E-23,	0.000000000000882096
                );
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
      if Tj[j] <> 0 then
        Result[i-1, j-1] := LiquidCp_a[1] * (LiquidCp_a[2] + LiquidCp_a[3] * LiquidCp_R[i]) * Tj[j] / (Tcc[i] + 273.15)
          + (LiquidCp_a[4] + LiquidCp_a[5] * LiquidCp_R[i]) * sqr(sqr(Tj[j] / (Tcc[i] + 273.15))) * Tj[j] / (Tcc[i] + 273.15)
          + LiquidCp_a[6] * sqr(LiquidCp_R[i]) / sqr(Tj[j] / (Tcc[i] + 273.15)) + LiquidCp_a[7] * LiquidCp_R[i]
          / (sqr(Tj[j] / (Tcc[i] + 273.15)) * Tj[j] / (Tcc[i] + 273.15))
          + LiquidCp_a[8] * sqr(sqr(Tj[j] / (Tcc[i] + 273.15))) * Tj[j] / (Tcc[i] + 273.15)
          + LiquidCp_k[i] * (LiquidCp_b[1] + LiquidCp_b[2] * sqr(Tj[j] / (Tcc[i] + 273.15)) +
          LiquidCp_b[3] * sqr(sqr(Tj[j] / (Tcc[i] + 273.15))) * Tj[j] / (Tcc[i] + 273.15))
          + sqr(LiquidCp_k[i]) * (LiquidCp_b[4] + LiquidCp_b[5] * sqr(Tj[j] / (Tcc[i] + 273.15))) + IdealGasCompCp[i-1, j-1]
      else
        Result[i-1, j-1] := 0;
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
  CompIntLiqCp: TArrOfArrOfDouble;
  CompIntIdGasCp: TArrOfArrOfDouble;
  LiquidCompHeatCapasity: TArrOfArrOfDouble;
  Tsat: arrComp;
  Tdp: arrComp;
  s: double;

begin
  SetLength(compH_v, NComp, NTrays);
  SetLength(compH_l, NComp, NTrays);
  SetLength(CompIntLiqCp, NComp, NTrays);
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
            CompIntLiqCp[i-1, j-1] := IntegralLiquidCp(Tsat, Tj[j])[i];
        end;
    end;

  for j := 1 to NTrays do
    begin
      s := 0;
      for i := 1 to NComp do
        s := s + {(dHf298[i] + (CompIntIdGasCp[i-1, j-1] - dHvap[i]
          + CompItnLiqCp[i-1, j-1]) * 4.1868)} LiquidCompHeatCapasity[i-1, j-1] * 4.1868 * (Tsat[i] - Tj[j]) * xij[i-1, j-1];
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
      if Tj[j] = 0 then
        continue;
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
  Qj: arrTrays; ej: arrTrays; LD, VD: double; rB: double; var Vj: arrTrays);
var
  j, k: integer;
  alp, bet, gam: arrTrays;
  s: double;

begin
  Vj[NTrays] := {Fj[NTrays] + Lj[NTrays-1] - Uj[NTrays] - Wj[NTrays]{Lj[NTrays-1] * 1.15}rB * Lj[NTrays];
  Vj[2] := Lj[1] + LD + VD;
  for j := NTrays-1 downto 2 do
    begin
      s := 0;
      for k := 1 to j-1 do
        s := s + (Fj[k] - Uj[k] - Wj[k]);
      gam[j] := (Hf_l[j] * (1 - ej[j]) + Hf_v[j] * ej[j] - H_l[j]) * Fj[j]
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

function TMatBalance.teta_method(Fj: arrTrays; zf, xij: TArrOfArrOfDouble; D: Double; B: Double; Dco: double): double;
const
  eps = 1e-5;
var
  s: double;
  i, j: double;
  teta, r_teta: double;
  _a, _b: double;

  function g(teta: double): double;
  var
    i, j: integer;
    s, s1: double;
    tet: arrComp;
  begin
    s := 0;
    for i := 1 to NComp do
      begin
        s1 := 0;
        for j := 1 to NTrays do
          begin
            s1 := s1 + Fj[j] * zf[i-1, j-1];
          end;
        if xij[i-1, 0] <> 0 then
          tet[i] := s1 / (1 + teta * (B * xij[i-1, NTrays-1]) / (D * xij[i-1, 0]))
        else
          tet[i] := 0;
        s := s + tet[i];
      end;
    Result := s - Dco;
  end;

  function g1(teta: double): double;
  var
    i, j: integer;
    s, s1: double;
    tet: arrComp;
  begin
    s := 0;
    for i := 1 to NComp do
      begin
        s1 := 0;
        for j := 1 to NTrays do
          begin
            s1 := s1 + Fj[j] * zf[i-1, j-1];
          end;
        tet[i] := (B * xij[i-1, NTrays-1]) / (D * xij[i-1, 0]) * s1
          / sqr(1 + teta * (B * xij[i-1, NTrays-1]) / (D * xij[i-1, 0]));
        s := s + tet[i];
      end;
    Result := -s;
  end;

  function get_borders(x: double): double;
  const
    h = 1;
  begin
    while g(x) * g(x + h) / abs(g(x + h)) >= 0 do
      x := x + h;
    Result := x + h;
  end;

begin
  r_teta := 0;
  //teta := r_teta - g(r_teta) / g1(r_teta);
  _a := 0;
  _b := get_borders(_a);
  if g(_a) * g(_b) / abs(g(_b)) < 0 then
    begin
      repeat
        {r_teta := teta;
        teta := r_teta - g(r_teta) / g1(r_teta); }
        teta := _a + (_b - _a) / 2;
        if g(_a) * g(teta) < 0 then
          _b := teta
        else
          _a := teta
      until (abs(_a - _b) <= eps) or (g(teta) = 0);
      teta := (_a + _b) / 2;
      Result := teta;
    end
  else
    ShowMessage('Teta method, no solutions!');
end;

procedure TMatBalance.CalculateSideFlows_and_LiquidFlows(dj, Fj: arrTrays; ej: arrTrays;
  Vj: arrTrays; Lj0: arrTrays; LD, L0: double; zf, xij: TArrOfArrOfDouble; var Wj: arrTrays; var Uj: arrTrays; var Lj: arrTrays);
var
  {dj: arrTrays;}
  Rj: arrTrays;
  j: integer;
  rD, rB: double;
  B, D: double;
  s, s1, sd, sb: double;
  qj: arrTrays;
  teta: double;
  di, bi: arrComp;
  i: Integer;

begin
  for j := 1 to NTrays do
    begin
      {if Lj0[j] + Vj[j] <> 0 then
        dj[j] := (Uj[j] + Wj[j]) / (Lj0[j] + Vj[j])
      else
        dj[j] := 0; }
      if j = 1 then
        Rj[j] := (dj[j] / (dj[j] + 1)) * (Fj[j] + Vj[j+1] + Lj0[j-1])
    else
        Rj[j] := (dj[j] / (dj[j] + 1)) * (Fj[j] + Vj[j+1] + Lj0[j-1]);
      {Wj[j] := ej[j] * Rj[j];
      Uj[j] := (1 - ej[j]) * Rj[j];}
    end;

  Uj[1] := 1313;
  Uj[NTrays] := Fj[FeedTray1] + Fj[FeedTray2] - Uj[1];

  Lj[1] := L0;
  for j := 2 to NTrays-1 do
    Lj[j] := Fj[j] + Vj[j+1] - Vj[j] + Lj[j-1] - Rj[j];
  Lj[NTrays] := Uj[NTrays];


  rD := Lj[1] / LD;
  rB := Vj[NTrays] / Uj[NTrays];

  for j := 1 to NTrays do
    qj[j] := 1; //�� ���� ���� ��� ������� ����� ��������� � ����������� �� ���������

  s := 0;
  for j := 2 to NTrays-1 do
    s := s + (qj[j] + rD) * Fj[j] + (rD + 1) * Uj[j] - Wj[j];
  B := (rD * Fj[1] + (rD + 1) * Fj[NTrays] + s) / (rD + rB + 1){9.3};

  s := 0;
  for j := 2 to NTrays-1 do
    s := s + (rB + 1 - qj[j]) * Fj[j] + (rD + 1) * Uj[j] + Wj[j];
  D := ((rB + 1) * Fj[1] + rD * Fj[NTrays] + s) / (rD + rB + 1){4.5};

  teta := teta_method(Fj, zf, xij, D, B, 1313);
  sd := 0;
  sb := 0;
  for i := 1 to NComp do
    begin
      s1 := 0;
      for j := 1 to NTrays do
        s1 := s1 + Fj[j] * zf[i-1, j-1];
      if xij[i-1, 0] <> 0 then
        begin
          di[i] := s1 / (1 + teta * (B * xij[i-1, NTrays-1]) / (D * xij[i-1, 0]));
          bi[i] := teta * (B * xij[i-1, NTrays-1]) / (D * xij[i-1, 0]) * di[i];
        end;
      sd := sd + di[i];
      sb := sb + bi[i];
    end;
  D := sd;
  B := sB;

  Uj[NTrays] := B;
  Uj[1] := D;
  Lj[1] := D * rD;
  Lj[Ntrays] := B;

  s := 0;
  for j := 2 to NTrays-1 do
    s := s + (qj[j] + rD) * Fj[j] + (rD + 1) * Uj[j] - Wj[j];
  rB := (rD + Fj[1] + (rD + 1) * Fj[NTrays] + s) / B - rD - 1;

  Vj[NTrays] := B * rB;
  
  for j := 2 to NTrays-1 do
    begin
      s := 0;
      for i := 1 to j do
        s := s + Fj[i] - Uj[i] - Wj[i];
      Lj[j] := Vj[j+1] + s - Vj[1];
    end;

end;

procedure TMatBalance.LV_correction(rD: Double; rB: Double; Fj: arrTrays; Uj: arrTrays;
  Wj: arrTrays; var Lj: arrTrays; var Vj: arrTrays);
var
  j, k: integer;
  s: double;
begin
  Vj[NTrays] := rB * Uj[NTrays];
  Vj[2] := (1 + rD) * Uj[1] + Wj[1];
  Lj[NTrays] := Uj[NTrays];
  Lj[NTrays-1] := Vj[NTrays] + Lj[NTrays] + Wj[NTrays];

  for j := 2 to NTrays-1 do
    begin
      s := 0;
      for k := 1 to j do
        s := s + (Fj[k] - Uj[k] - Wj[k]);
      Lj[j] := Vj[j+1] + s - Vj[1];
    end;
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

procedure TMatBalance.correction_with_tray_efficiencies(tray_efficiency: arrTrays; Kij: TArrOfArrOfDouble;
  var xij: TArrOfArrOfDouble; var yij: TArrOfArrOfDouble);
var
  i, j: integer;
begin
  for j := 1 to NTrays-1 do
    for i := 0 to NComp-1 do
      xij[i, j] := xij[i, j-1] - tray_efficiency[j] * (xij[i, j-1] - xij[i, j]);
  for j := 0 to NTrays-2 do
    for i := 0 to NComp-1 do
      yij[i, j] := yij[i, j+1] + tray_efficiency[j] * (yij[i, j] - yij[i, j+1]);

  for i := 0 to Ncomp-1 do
    begin
      xij[i, 0] := yij[i, 0] / Kij[i, 0];
      yij[i, NTrays-1] := xij[i, NTrays-1] * Kij[i, NTrays-1];
    end;

end;

procedure TMatBalance.MatBalCalculation(Fl: arrTrays; Fv: arrTrays; Wl: arrTrays;
  Wv: arrTrays; T1: Double; TN: Double; P1: Double; PN: Double;
  D: double; LD: double;
  var Tj: arrTrays; var Lj: arrTrays; var Vj: arrTrays;
  var xij: TArrOfArrOfDouble; var yij: TArrOfArrOfDouble;
  var calcTj: TArrOfArrOfDouble;
  var calcLj: TArrOfArrOfDouble;
  var calcVj: TArrOfArrOfDouble;
  var n: integer);
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
  {xij: TArrOfArrOfDouble;}
  norm_xij: TArrOfArrOfDouble;
 { Tj: arrTrays; }
  Tfj: arrTrays;
  Pfj: arrTrays;
 { yij: TArrOfArrOfDouble;}
  H_l, H_v: arrTrays;
  xf, yf: TArrOfArrOfDouble;
  ej: arrTrays;
  ej_tr: arrTrays;
  Hf_l, Hf_v: arrTrays;
  {Lj: arrTrays; }
  {Vj: arrTrays; }
  Qj: arrTrays;
  Qc, Qr: Double;
  {calcTj: TArrOfArrOfDouble;}
  {calcLj: TArrOfArrOfDouble;
  calcVj: TArrOfArrOfDouble;
  n: integer; }
  omegaj: ArrTrays; // ���� ���� � ������� �������
  TrayMaterialBalanceError: arrTrays;
  TrayHeatBalanceError: arrTrays;
  dj: arrTrays;
  rB: double;
  tray_efficiency: arrTrays;

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
  {
  zf[0, FeedTray1-1] := 2.22222222222222e-002;
  zf[1, FeedTray1-1] := 4.44444444444444e-002;
  zf[2, FeedTray1-1] := 6.66666666666667e-002;
  zf[3, FeedTray1-1] := 8.88888888888889e-002;
  zf[4, FeedTray1-1] := 0.111111111111111;
  zf[5, FeedTray1-1] := 0.133333333333333;
  zf[6, FeedTray1-1] := 0.155555555555556;
  zf[7, FeedTray1-1] := 0.177777777777778;
  zf[8, FeedTray1-1] := 0.200000000000000;
  }
  zf[ 0, FeedTray1-1] := 0;
  zf[ 1, FeedTray1-1] := 0;
  zf[ 2, FeedTray1-1] := 0;
  zf[ 3, FeedTray1-1] := 0.532434510310423;
  zf[ 4, FeedTray1-1] := 0.139469584691842;
  zf[ 5, FeedTray1-1] := 0.0964475485590856;
  zf[ 6, FeedTray1-1] := 0.0509975946839561;
  zf[ 7, FeedTray1-1] := 0.0430189661776558;
  zf[ 8, FeedTray1-1] := 0.011361250705204;
  zf[ 9, FeedTray1-1] := 0.00924679155267898;
  zf[10, FeedTray1-1] := 0.00799659213845727;
  zf[11, FeedTray1-1] := 0;
  zf[12, FeedTray1-1] := 0.00894327768172518;
  zf[13, FeedTray1-1] := 0.000547052895275276;
  zf[14, FeedTray1-1] := 0.000530266973209197;
  zf[15, FeedTray1-1] := 0.00450934263732591;
  zf[16, FeedTray1-1] := 0.00037493268704121;
  zf[17, FeedTray1-1] := 0;
  zf[18, FeedTray1-1] := 0.00723621406102272;
  zf[19, FeedTray1-1] := 0;
  zf[20, FeedTray1-1] := 0.0109478200055871;
  zf[21, FeedTray1-1] := 0.0195150395003681;
  zf[22, FeedTray1-1] := 0.0249593124719107;
  zf[23, FeedTray1-1] := 0.012700891777383;
  zf[24, FeedTray1-1] := 0;
  zf[25, FeedTray1-1] := 0.00295291319519114;
  zf[26, FeedTray1-1] := 0.00138658785404103;
  zf[27, FeedTray1-1] := 0.00054943739468231;
  zf[28, FeedTray1-1] := 0.00336017803832046;
  zf[29, FeedTray1-1] := 0.0105138940076147;
  zf[30, FeedTray1-1] := 0;

  zf[ 0, FeedTray2-1] := 0.00075792344443219;
  zf[ 1, FeedTray2-1] := 0;
  zf[ 2, FeedTray2-1] := 0.00218269096480515;
  zf[ 3, FeedTray2-1] := 0.95590636432036;
  zf[ 4, FeedTray2-1] := 0.0362677562997908;
  zf[ 5, FeedTray2-1] := 0;
  zf[ 6, FeedTray2-1] := 0;
  zf[ 7, FeedTray2-1] := 0;
  zf[ 8, FeedTray2-1] := 0;
  zf[ 9, FeedTray2-1] := 0;
  zf[10, FeedTray2-1] := 0;
  zf[11, FeedTray2-1] := 0;
  zf[12, FeedTray2-1] := 0;
  zf[13, FeedTray2-1] := 0;
  zf[14, FeedTray2-1] := 0;
  zf[15, FeedTray2-1] := 0;
  zf[16, FeedTray2-1] := 0;
  zf[17, FeedTray2-1] := 0;
  zf[18, FeedTray2-1] := 0;
  zf[19, FeedTray2-1] := 0;
  zf[20, FeedTray2-1] := 0;
  zf[21, FeedTray2-1] := 0.00488526497061137;
  zf[22, FeedTray2-1] := 0;
  zf[23, FeedTray2-1] := 0;
  zf[24, FeedTray2-1] := 0;
  zf[25, FeedTray2-1] := 0;
  zf[26, FeedTray2-1] := 0;
  zf[27, FeedTray2-1] := 0;
  zf[28, FeedTray2-1] := 0;
  zf[29, FeedTray2-1] := 0;
  zf[30, FeedTray2-1] := 0;

  for j := 1 to NTrays do
    begin
      Fj[j] := 0;
      Wj[j] := 0;
      Uj[j] := 0;
      Qj[j] := 0;
      tray_efficiency[j] := 1;
    end;
  Fj[FeedTray1] := {13.8}1618;
  Fj[FeedTRay2] := 103.1;
  Pj[1] := P1;
  Pj[Ntrays] := PN;
  for j := 2 to NTrays-1 do
    Pj[j] := Pj[j-1] + (Pj[NTrays] - Pj[1]) / NTrays;

  for j := 1 to NTrays do
    begin
      Tfj[j] := 0;
      Pfj[j] := Pj[j];
      omegaj[j] := 0; // �.�. ��� ����� � ���� ��������
    end;
  Tfj[FeedTray1] := 40 + 273.15;
  Tfj[FeedTray2] := 45 + 273.15;
  Pfj[FeedTray1] := 0.65;
  Pfj[FeedTray2] := 0.65;

  Tj_0 := getTj_0(Fj, zf);
  InitGuessLV(Fj, 0, 1313, 376.2, Lj0, Vj0, Uj, Wj, dj, rB);
  WilsonCorrelation(Tcc, Pcc, omega, NTrays, Tj_0, Pj, Kij);
  CalculateLiquidMoleFractions(Fj, Lj0, Vj0, Uj, Wj, zf, Kij, xij);
  xij := Normalize(xij);
  Secant(50, 900, Tj_0, Pj, xij, yij, Tj);
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
    Secant(50, 500, Tj_0, Pj, xij, yij, Tj);
    WilsonCorrelation(Tcc, Pcc, omega, NTrays, Tj, Pj, Kij);
    CalculateVaporMoleFractions(xij, Kij, yij);
    yij := Normalize(yij);
    correction_with_tray_efficiencies(tray_efficiency, Kij, xij, yij);
    xij := Normalize(xij);
    yij := Normalize(yij);
    CalculateEnthalpies(Tj, Pj, xij, yij, H_l, H_v);
    getEquilibrium(zf, Tfj, Pj, xf, yf, ej);
    CalculateEnthalpies(Tfj, Pfj, xf, yf, Hf_l, Hf_v);
    CalculateVaporFlowRates(Hf_l, Hf_v, H_l, H_v, Fj, Lj0, Uj, Wj, Qj, ej, 1313, 0, rB, Vj);
    CalculateHeatDuties(Fj, ej, Lj0, Vj, Uj, Wj, Hf_l, Hf_v, H_l, H_v, Qc, Qr);
    ej_tr := getTrays_ej(Fj, Lj, Vj, zf, xij, yij, Tj, Pj); // ��� ���������, ��� �� �����
    CalculateSideFlows_and_LiquidFlows(dj, Fj, omegaj, Vj, Lj0, 1313, 376.2, zf, xij, Wj, Uj, Lj);
    LV_correction(376.2 / 1313, rB, Fj, Uj, Wj, Lj, Vj);
    TrayMaterialBalanceError := getTrayMaterialBalanceError(Fj, Uj, Wj, Lj, Vj);
    TrayHeatBalanceError := getTrayHeatBalanceError(Fj, Uj, Wj, Lj, Vj, ej, Hf_l, Hf_v, H_l, H_v);
    n := n + 1;
    if n >= 1e3 then
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
        Tj_0[j] := Tj[j];
      end;

  until {getErrorValue(n, calcTj, calcLj, calcVj)}(abs(Uj[1] - D)) + (abs(Lj[1] - LD)) <= tolerance;
  ShowMessage(IntToStr(n))


  end;

end.
