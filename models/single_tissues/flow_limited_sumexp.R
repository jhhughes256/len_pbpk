# Flow limited model for use with single tissues
# -----------------------------------------------------------------------------
# Default tissue values are for heart
# Uses cubic forcing function to feed in arterial concentrations
# Default values for sum of exponentials are not particularly useful

  code <- '
$INIT
  Cven = 0  // Venous Blood

$PARAM
  // Regional blood flow (mL/min)
  Q = 0.9227,

  // Apparent volume of distribution (mL)
  V = 0.1,

  // Real volume of distribution (mL)
  Vreal = 0.125,

  // Arterial forcing function (sum of exponentials)
  M1 = -0.017,
  M2 = -2.2,
  B = 5.3

$MAIN
  double Vrat = V/Vreal;
  double Cart = exp(M1*TIME + B) - exp(M2*TIME + B);

$ODE
  dxdt_Cven = Q*(Cart - Cven)/V;

$TABLE
  double Ctis = Cven*Vrat;

$CAPTURE
  Cart Cven Ctis Q V M1 M2 B
  '
# Compile the model code
  flowmod <- mcode("flowlim", code)
