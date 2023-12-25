unit TypesUnit;

interface

type
  TArrOfDouble = array of Double;
  TMatrixOfDouble = array of array of Double;
  TFoo = function(x: TArrOfDouble): Double;
  TCubicEquationMethod = function(a, b, c: Double): TArrOfDouble;
  TObjectiveFunction = reference to function(x: Double): Double;

implementation

end.
