# Membrane limited model for use with single tissues
# -----------------------------------------------------------------------------
# Default tissue values are for heart
# Uses cubic forcing function to feed in arterial concentrations
# Default values for linear equation are not particularly useful

code <- '
$INIT
  Cart = 0,  // Arterial blood
  Cven = 0,  // Venous Blood
  Cmem = 0  // Tissue

$PARAM
  // Regional blood flow (mL/min)
  Q = 0.9227,

  // Apparent volume of distribution (mL)
  V1 = 0.1,
  V2 = 0.15,

  // Permeability
  PS = 0.01,

  // Real volume of distribution (mL)
  Vreal = 0.125,

  // Arterial forcing function (linear)
  M = 0

$MAIN
  double Vrat = (V1 + V2)/Vreal;

$ODE
  dxdt_Cart = M;
  dxdt_Cven = (Q*(Cart - Cven) + PS*(Cmem - Cven))/V1;
  dxdt_Cmem = PS*(Cven - Cmem)/V2;

$TABLE
  double Ctis = Cmem*Vrat;

$CAPTURE
  Cart Cven Cmem Ctis Q V1 PS V2 M
'
# Compile the model code
  membmod <- mcode("memblim", code)
