# Mouse model for lenalidomide
# -----------------------------------------------------------------------------
# Main source for reference values:
#   Brown RP, Delp MD, Lindstedt SL, Rhomberg LR, Beliles RP. Physiological
#   Parameter Values for Physiologically Based Pharmacokinetic Models.
#   Toxicology and Industrial Health. 1997;13(4):407-84.
# Source used for brain blood flow:
#   Davies B, Morris T. Physiological Parameters in Laboratory Animals and
#   Humans. Pharmaceutical Research. 1993;10(7):1093-5.
# -----------------------------------------------------------------------------
# Model specifications
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Build: Feed spleen into liver
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Model code
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  code <- '
$INIT
  // Initial conditions for compartments
  Apa = 0,  // Vascular mixing
  Aart = 0,  // Arterial Blood
  Alng = 0,  // Lung
  Ahrt = 0,  // Heart
  Akid = 0,  // Kidney
  Alvr = 0,  // Liver
  Abra = 0,  // Brain
  Amsc = 0,  // Muscle
  Aspl = 0,  // Spleen sinus
  Abod = 0  // Rest of body

$PARAM
  // Standard Physiological Parameters (Brown et. al 1997)
  WTstd = 25,  // Weight (g)
  COstd = 13.98,  // Cardiac Output (ml/min)

  // Regional Blood Flow (mL/min)
  Qhrtstd = 0.9227,
  Qkidstd = 1.272,
  Qlvrstd = 2.251,
  Qbrastd = 0.4613,
  Qmscstd = 2.223,
  Qsplstd = 0.1580,  // Imputed from (Davies et. al 1993)
  PSlngstd = 0.2,
  PSsplstd = 0.001,

  // Tissue Mass Balance (mL)
  Vmixstd = 1.225,
  Vlngstd = 0.175,
  Vhrtstd = 0.125,
  Vkidstd = 0.425,
  Vlvrstd = 1.375,
  Vbrastd = 0.425,
  Vmscstd = 9.6,
  Vsplstd = 0.0875

  // Individual Covariate Values
  WT = 20

$MAIN
  // Remainder of cardiac output and volume
  double Qbodstd = COstd-(Qhrtstd+Qkidstd+Qlvrstd+Qbrastd+Qmscstd);
  double Vbodstd = WTstd-(Vbrastd+Vlvrstd+Vmscstd+Vhrtstd+Vsplstd+Vkidstd);

  // Allometric scaling of blood flows, clearances and permeabilities
  double PSlng = PSlngstd*pow(WT/WTstd,0.75);
  double Qhrt = Qhrtstd*pow(WT/WTstd,0.75);
  double Qkid = Qkidstd*pow(WT/WTstd,0.75);
  double Qlvr = Qlvrstd*pow(WT/WTstd,0.75);
  double Qbra = Qbrastd*pow(WT/WTstd,0.75);
  double Qmsc = Qmscstd*pow(WT/WTstd,0.75);
  double Qspl = Qsplstd*pow(WT/WTstd,0.75);
  double Qbod = Qbodstd*pow(WT/WTstd,0.75);
  double Qco = Qhrt+Qkid+Qlvr+Qbra+Qmsc+Qbod;

  // Apparent distribution volumes with allometric scaling
  double Vmix = Vmixstd*pow(WT/WTstd,1);
  double Vlng = Vlngstd*pow(WT/WTstd,1);
  double Vhrt = Vhrtstd*pow(WT/WTstd,1);
  double Vkid = Vkidstd*pow(WT/WTstd,1);
  double Vlvr = Vlvrstd*pow(WT/WTstd,1);
  double Vbra = Vbrastd*pow(WT/WTstd,1);
  double Vmsc = Vmscstd*pow(WT/WTstd,1);
  double Vspl = Vsplstd*pow(WT/WTstd,1);
  double Vbod = Vbodstd*pow(WT/WTstd,1);

$ODE
  dxdt_Apa = -Qco*Apa/Vmix +Qhrt*Ahrt/Vhrt +Qkid*Akid/Vkid +Qlvr*Alvr/Vlvr
    +Qbra*Abra/Vbra +Qmsc*Amsc/Vmsc +Qbod*Abod/Vbod;
  dxdt_Aart = Qco*(Apa/Vmix -Aart/Vlng) +PSlng*(Alng -Aart);
  dxdt_Alng = PSlng*(Aart -Alng);
  dxdt_Ahrt = Qhrt*(Aart/Vlng -Ahrt/Vhrt);
  dxdt_Akid = Qkid*(Aart/Vlng -Akid/Vkid);
  dxdt_Alvr = (Qlvr-Qspl)*Aart/Vlng -Qlvr*Alvr/Vlvr +Qspl*Aspl/Vspl;
  dxdt_Abra = Qbra*(Aart/Vlng -Abra/Vbra);
  dxdt_Amsc = Qmsc*(Aart/Vlng -Amsc/Vmsc);
  dxdt_Aspl = Qspl*(Aart/Vlng -Aspl/Vspl);
  dxdt_Abod = Qbod*(Aart/Vlng -Abod/Vbod);

$TABLE  // Determine individual predictions
  double Cpa = Apa/Vmix;
  double Cart = Aart/Vlng;
  double Clng = Alng/Vlng;
  double Chrt = Ahrt/Vhrt;
  double Ckid = Akid/Vkid;
  double Clvr = Alvr/Vlvr;
  double Cbra = Abra/Vbra;
  double Cmsc = Amsc/Vmsc;
  double Cspl = Aspl/Vspl;
  double Cbod = Abod/Vbod;

$CAPTURE
  Cpa Cart Clng Chrt Ckid Clvr Cbra Cmsc Cspl Cbod
  COstd WTstd Qhrt Qkid Qlvr Qbra Qmsc Qspl Qbod Qco
  Vmix Vlng Vhrt Vkid Vlvr Vbra Vmsc Vspl Vbod
'
# Compile code
  brown.mod <- mcode("mouseBROWN", code)
