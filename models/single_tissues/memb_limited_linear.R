# Membrane limited model for use with single tissues
# -----------------------------------------------------------------------------
# Default tissue values are for heart
# Uses cubic forcing function to feed in arterial concentrations
# Default values for cubic spline are not particularly useful

code <- '
$INIT
  Cven = 0,  // Venous Blood
  Ctis = 0  // Tissue

$PARAM
  // Regional blood flow (mL/min)
  Q = 0.9227,

  // Apparent volume of distribution (mL)
  V1 = 0.125,
  V2 = 0.1,

  // Permeability
  PS = 0.01,

  // Real volume of distribution (mL)
  Vreal = 0.1,

  // Arterial forcing function (linear)
  M = -0.1,
  B = 7

$MAIN
  // double Vrat = V/Vreal;
  double Cart = M*TIME + B;

$ODE
  dxdt_Cven = (Q*(Cart - Cven) + PS*(Ctis - Cven))/V1;
  dxdt_Ctis = PS*(Cven - Ctis)/V2;

$TABLE
  // double Ctis = Cven*Vrat;

$CAPTURE
  Cart Cven Ctis Q V1 PS V2 M B
'
# Compile the model code
  membmod <- mcode("memblim", code)
