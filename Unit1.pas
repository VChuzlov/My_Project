unit Unit1;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls, UMaterialBalance, UPhase_Equilibrium,
  Vcl.OleServer, ExcelXP;

type
  TForm1 = class(TForm)
    Button1: TButton;
    EA1: TExcelApplication;
    Button2: TButton;
    procedure get_report(Lj, Vj, Tj: arrTrays; xij, yij: TArrOfArrOfDouble;
                         calcTj: TArrOfArrOfDouble;
                         calcLj: TArrOfArrOfDouble;
                         calcVj: TArrOfArrOfDouble;
                         n: integer);
    procedure Button1Click(Sender: TObject);
    procedure Button2Click(Sender: TObject);
    procedure FormClose(Sender: TObject; var Action: TCloseAction);


  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  Form1: TForm1;
  Fl, Fv, Uj, Wj, Lj, Vj, Tj: arrTrays;
  yij, xij: TArrOfArrOfDouble;
  calcTj: TArrOfArrOfDouble;
  calcLj: TArrOfArrOfDouble;
  calcVj: TArrOfArrOfDouble;
  n: integer;

implementation

{$R *.dfm}

procedure TForm1.get_report(Lj: arrTrays; Vj: arrTrays; Tj: arrTrays;
                            xij: TArrOfArrOfDouble; yij: TArrOfArrOfDouble;
                            calcTj: TArrOfArrOfDouble;
                            calcLj: TArrOfArrOfDouble;
                            calcVj: TArrOfArrOfDouble;
                            n: integer);
var
  i, j: integer;
begin
  EA1.Connect;
  EA1.Workbooks.Add(Null, 0);
  EA1.Cells.NumberFormat:= '@';
  EA1.Cells.Item[1, 1] := 'Tj';
  EA1.Cells.Item[1, 2] := 'Lj';
  EA1.Cells.Item[1, 3] := 'Vj';
  for j := 1 to NTrays do
    begin
      EA1.Cells.Item[j+1, 1]:= Tj[j];
      EA1.Cells.Item[j+1, 2]:= Lj[j];
      EA1.Cells.Item[j+1, 3]:= Vj[j];
    end;
  for i := 1 to NComp do
    for j := 1 to NTrays do
      begin
        EA1.Cells.Item[i+NTrays+2, j] := xij[i-1, j-1];
        EA1.Cells.Item[i+2*NTrays, j] := yij[i-1, j-1];
      end;
  for i := 1 to n do
    for j := 1 to NTrays do
      begin
        EA1.Cells.Item[i+3*NTrays, j] := calcTj[i-1, j-1];
        EA1.Cells.Item[i+3*NTrays+n+1, j] := calcLj[i-1, j-1];
        EA1.Cells.Item[i+3*NTrays+2*n+2, j] := calcVj[i-1, j-1];
      end;
  EA1.Visible[0]:= True;
end;

procedure TForm1.Button2Click(Sender: TObject);
begin
  get_report(Lj, Vj, Tj, xij, yij, calcTj, calcLj, calcVj, n)
end;

procedure TForm1.Button1Click(Sender: TObject);
var
  MatBal: TMatBalance;

begin
  MatBal.Create;
  MatBal.MatBalCalculation(Fl, Fv, Uj, Wj, 48, 120, 0.576, 0.655, 1313, 376.2, Tj, Lj, Vj,
                           xij, yij, calcTj, calcLj, calcVj, n);
  ShowMessage('Success !!!')
end;

procedure TForm1.FormClose(Sender: TObject; var Action: TCloseAction);
begin
  EA1.Disconnect;
  EA1.Free
end;

end.
