# Flow limited model for use with single tissues
# -----------------------------------------------------------------------------
# Default tissue values are for heart
# Uses cubic forcing function to feed in arterial concentrations
# Default values for cubic spline are not particularly useful

  code <- '
$INIT
  Cven = 0  // Venous Blood

$PARAM
  // Regional blood flow (mL/min)
  Q = 0.9227,

  // Apparent volume of distribution (mL)
  V = 0.125,

  // Real volume of distribution (mL)
  Vreal = 0.1,

  // Arterial forcing function (linear)
  M = -0.1,
  B = 7

$MAIN
  double Vrat = V/Vreal;
  double Cart = M*TIME + B;

$ODE
  dxdt_Cven = Q*(Cart - Cven)/V;

$TABLE
  double Ctis = Cven*Vrat;

$CAPTURE
  Cart Cven Ctis Q V  M B
  '
# Compile the model code
  flowmod <- mcode("flowlim", code)
