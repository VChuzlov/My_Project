program Rectification;

uses
  Vcl.Forms,
  Unit1 in 'Unit1.pas' {Form1},
  UMaterialBalance in 'UMaterialBalance.pas',
  UPhase_Equilibrium in 'UPhase_Equilibrium.pas';

{$R *.res}

begin
  Application.Initialize;
  Application.MainFormOnTaskbar := True;
  Application.CreateForm(TForm1, Form1);
  Application.Run;
end.
