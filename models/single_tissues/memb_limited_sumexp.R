# Membrane limited model for use with single tissues
# -----------------------------------------------------------------------------
# Default tissue values are for heart
# Uses cubic forcing function to feed in arterial concentrations
# Default values for sum of exponentials are not particularly useful

code <- '
$INIT
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

  // Arterial forcing function (sum of exponentials)
  M1 = -0.017,
  M2 = -2.2,
  B = 5.3

$MAIN
  double Vrat = (V1 + V2)/Vreal;
  double Cart = exp(M1*TIME + B) - exp(M2*TIME + B);

$ODE
  dxdt_Cven = (Q*(Cart - Cven) + PS*(Cmem - Cven))/V1;
  dxdt_Cmem = PS*(Cven - Cmem)/V2;

$TABLE
  double Ctis = Cmem*Vrat;

$CAPTURE
  Cart Cven Cmem Ctis Q V1 PS V2 M1 M2 B
'
# Compile the model code
  membmod <- mcode("memblim", code)
