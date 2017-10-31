# Flow limited model for use with single tissues
# -----------------------------------------------------------------------------
# Default tissue values are for heart
# Uses cubic forcing function to feed in arterial concentrations
# Default values for linear equation are not particularly useful

  code <- '
$INIT
  Cart = 0  // Arterial Blood

$PARAM
  // Arterial forcing function (linear)
  M = 0

$ODE
  dxdt_Cart = M;

$CAPTURE
  Cart M
  '
# Compile the model code
  linmod <- mcode("linear", code)
