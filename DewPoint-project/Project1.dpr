program Project1;

uses
  Vcl.Forms,
  Unit1 in 'Unit1.pas' {Form1},
  DewPointUnit in 'DewPointUnit.pas',
  FunctionsUnit in 'FunctionsUnit.pas',
  TypesUnit in 'TypesUnit.pas';

{$R *.res}

begin
  Application.Initialize;
  Application.MainFormOnTaskbar := True;
  Application.CreateForm(TForm1, Form1);
  Application.Run;
end.
