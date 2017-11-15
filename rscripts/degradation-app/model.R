# Degradation example model
# -----------------------------------------------------------------------------
# A proposed model for drug degradation and movement for protein binding 
# experiment using a vertical dialysis apparatus
# Units are:
# - drug amount = nanomoles
# - volume = millilitres (mL)
# - conc = nanomoles per millilitre (nM) or micromoles per litre (uM)
# - time = minutes (min)

  code <- '
$INIT
  Cpu = 0,  // Plasma Unbound
  Cpb = 0,  // Plasma Bound
  Cpd = 0,  // Plasma Degraded
  Cdu = 0,  // Dialysate Unbound
  Cdd = 0  // Dialysate Degraded

$PARAM
  // Dialysis membrane apparatus parameters
  Qm = 0.001,  // Membrane flow rate (mL/min)
  V = 0.1,  // Well volume (mL)

  // Plasma protein binding parameters
  fu = 0.7,  // Fraction unbound
  Cprot = 4500,
  Koff = 0.5,

  // Degradation rates (mL/min)
  CLp = 2*8.664*10**-3,  // plasma
  CLd = 2.888*10**-3  // dialysate

$MAIN
  double Kon = Koff*(1/fu-1)/Cprot;

$ODE
  dxdt_Cpu = Qm*(Cdu - Cpu)/V -Kon*Cprot*Cpu +Koff*Cpb -CLp*Cpu;
  dxdt_Cpb = Kon*Cprot*Cpu -Koff*Cpb;
  dxdt_Cpd = CLp*Cpu;
  dxdt_Cdu = Qm*(Cpu - Cdu)/V -CLd*Cdu;
  dxdt_Cdd = CLd*Cdu;

$CAPTURE
  Cpu Cpb Cpd Cdu Cdd 
'
# Compile the model code
  degmod <- mcode("lena-degradation", code)
