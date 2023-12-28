﻿unit Unit1;

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
    -82.45099487,  32.27800903, 96.74801025, 134.9460083,  152.04900513,
    187.24801025, 196.4500061, 234.74801025, 267.00802002, 295.44802246
  ];
  Pc: TArrOfDouble = [
    4640.68017578, 4883.85009766, 4256.66015625, 3647.62011719, 3796.62011719,
    3333.59008789, 3375.12011719, 3031.62011719, 2736.7800293,  2496.62011719
  ];
  Af: TArrOfDouble = [
    0.0114984,  0.0986,     0.1524,     0.18479,   0.20100001, 0.22224,
    0.25389001, 0.30070001, 0.34979001, 0.40180001
  ];
  Vc: TArrOfDouble = [
    0.0989999,  0.148,      0.2,        0.26300001, 0.25499001, 0.30799001,
    0.31099001, 0.368,      0.42598,    0.486
  ];
  Pressure: Double = 101.325;
  Yi: TArrOfDouble = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1];

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
