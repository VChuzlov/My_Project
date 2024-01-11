unit Unit1;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants,
  System.Classes, System.DateUtils,
  Vcl.Graphics, Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls,
  ConvertersUnit, TypesUnit, DewPointUnit;

type
  TForm1 = class(TForm)
    Button1: TButton;
    procedure Button1Click(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  Form1: TForm1;

implementation

{$R *.dfm}

procedure TForm1.Button1Click(Sender: TObject);
const
  Tc: TArrOfDouble = [
     -82.451,   32.278,   96.748,  134.946,  152.049,
     187.248,  196.450,  234.748,  267.008,  295.448,
     321.448,  344.448,  365.149,  385.149,  402.649,
     420.850,  433.850,  443.850,  460.220,  472.110,
     482.779,  494.850,  504.850,  513.850,  522.850,
     530.850,  538.850,  545.850,  552.850,  558.850
  ];
  Pc: TArrOfDouble = [
    4640.680, 4883.850, 4256.660, 3647.620, 3796.620,
    3333.590, 3375.120, 3031.620, 2736.780, 2496.620,
    2300.070, 2107.550, 1964.930, 1829.920, 1723.530,
    1620.180, 1516.810, 1420.560, 1316.900, 1213.470,
    1116.950, 1160.000, 1110.000, 1060.000, 1020.000,
     980.000,  950.000,  910.000,  883.000,  850.000
  ];
  Af: TArrOfDouble = [
          0.011,    0.099,    0.152,    0.185,    0.201,
         0.222,    0.254,    0.301,    0.350,    0.402,
         0.445,    0.488,    0.535,    0.562,    0.623,
         0.679,    0.706,    0.765,    0.770,    0.800,
         0.827,    0.907,    0.942,    0.972,    1.026,
         1.071,    1.105,    1.154,    1.214,    1.238
  ];
  Vc: TArrOfDouble = [
          0.099,    0.148,    0.200,    0.263,    0.255,
         0.308,    0.311,    0.368,    0.426,    0.486,
         0.543,    0.602,    0.660,    0.713,    0.780,
         0.830,    0.880,    0.940,    1.006,    1.070,
         1.130,    1.190,    1.245,    1.301,    1.356,
         1.411,    1.465,    1.518,    1.571,    1.622
  ];
  Pressure: Double = 101.325;
//  Yi: TArrOfDouble = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1];
    Yi: TArrOfDouble = [
      0.0333, 0.0333, 0.0333, 0.0333, 0.0333, 0.0333, 0.0333, 0.0333, 0.0333,
      0.0333, 0.0333, 0.0333, 0.0333, 0.0333, 0.0333, 0.0333, 0.0333, 0.0333,
      0.0333, 0.0333, 0.0333, 0.0333, 0.0333, 0.0333, 0.0333, 0.0333, 0.0333,
      0.0333, 0.0333, 0.0333
    ];
var
  uc: TUnitsConverter;
  dp: TDewPoint;
  TDew: Double;
  Start, Stop: TDateTime;
  Elapsed: Double;
  iCounterPerSec: TLargeInteger;
  T1, T2: TLargeInteger;  // значение счётчика ДО и ПОСЛЕ операции
begin
  Start := Now();
  QueryPerformanceFrequency(iCounterPerSec);  // определили частоту счётчика
  QueryPerformanceCounter(T1);  // засекли время начала операции

  uc := TUnitsConverter.Create();
  dp := TDewPoint.Create(
    Pressure, Yi, Tc, Pc, Af, Vc
  );
  TDew := dp.Calculation();

  QueryPerformanceCounter(T2);  // засекли время окончания
  Stop := Now();
  Elapsed := SecondsBetween(Start, Stop);  // время в секундах
  ShowMessage(
    FloatToStr(uc.Temperature.RankineToCelsius(TDew))
    + ' Elapsed at: ' + FormatFloat('0.0000', (T2 - T1)/iCounterPerSec) + ' sec.'
  );
end;

end.
