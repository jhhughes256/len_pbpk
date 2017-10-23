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
  PS = 0.1,

  // Real volume of distribution (mL)
  Vreal = 0.1,

  // Arterial forcing function (cubic)
  CUBT = 0,
  COF1 = 3,
  COF2 = 0.05,
  COF3 = -0.0042,
  COF4 = 0.000047

$MAIN
  // double Vrat = V/Vreal;
  double Tcub = TIME - CUBT;

$ODE
  double Cart = COF1 + COF2*Tcub + COF3*pow(Tcub,2) + COF4*pow(Tcub,3);
  dxdt_Cven = (Q*(Cart - Cven) + PS*(Ctis - Cven))/V1;
  dxdt_Ctis = PS*(Cven - Ctis)/V2;

$TABLE
  // double Ctis = Cven*Vrat;

$CAPTURE
  Cart Cven Ctis Q V1 PS V2 CUBT COF1 COF2 COF3 COF4
'
# Compile the model code
  membmod <- mcode("memblim", code)
