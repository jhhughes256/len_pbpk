# Mouse model for lenalidomide
# -----------------------------------------------------------------------------
# Main source for reference values:
#   Davies B, Morris T. Physiological Parameters in Laboratory Animals and
#   Humans. Pharmaceutical Research. 1993;10(7):1093-5.
# Source used for brain blood flow:
#   Brown RP, Delp MD, Lindstedt SL, Rhomberg LR, Beliles RP. Physiological
#   Parameter Values for Physiologically Based Pharmacokinetic Models.
#   Toxicology and Industrial Health. 1997;13(4):407-84.
# -----------------------------------------------------------------------------
# Model specifications
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Build: Adding GIT compartment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Model code
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  code <- '
$INIT
  // Initial conditions for compartments
  Apa = 0,  // Vascular mixing
  Aart = 0,  // Arterial Blood
  Alng = 0,  // Lung
  Abra = 0,  // Brain
  Alvr = 0,  // Liver
  Agut = 0,  // GIT
  Aspr = 0,  // Spleen sinus
  Asps = 0,  // Spleen pulp
  Akid = 0,  // Kidney
  Ahrt = 0,  // Heart
  Amsc = 0,  // Muscle
  Abod = 0  // Rest of body

$PARAM
  // Standard Physiological Parameters (Brown et. al 1997)
  WTstd = 20,  // Weight (g)
  COstd = 8.0,  // Cardiac Output (ml/min)
  GFRstd = 0.28,  // Glomerular Filtration Rate (ml/min)

  // Regional Blood Flow (mL/min)
  Qbrastd = 0.264,  // Imputed from Brown et. al 1997
  Qlvrstd = 1.8,
  Qgutstd = 1.5,
  Qsplstd = 0.09,
  Qkidstd = 1.3,
  Qhrtstd = 0.28,
  Qmscstd = 0.91,

  // Permeability Surface Area (mL/min)
  PSlngstd = 0.2,
  PSsplstd = 0.0158,

  // Tissue Mass Balance (mL)
  Vmixstd = 1.7,
  Vlngstd = 0.1,
  Vbrastd = 0.34,  // Imputed from Brown et. al 1997
  Vlvrstd = 1.3,
  Vgutstd = 1.5,
  Vsplstd = 0.1,
  Vkidstd = 0.34,
  Vhrtstd = 0.095,
  Vmscstd = 10,

  // Individual Covariate Values
  WT = 20

$MAIN
  // Remainder of cardiac output and volume
  double Qbodstd = COstd-(Qbrastd+Qlvrstd+Qkidstd+Qhrtstd+Qmscstd);
  double Vbodstd = WTstd-(Vbrastd+Vlvrstd+Vmscstd+Vhrtstd+Vsplstd+Vkidstd);

  // Allometric scaling of blood flows, clearances and permeabilities
  double PSlng = PSlngstd*pow(WT/WTstd,0.75);
  double PSspl = PSsplstd*pow(WT/WTstd,0.75);
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
  dxdt_Alvr = (Qlvr-Qspl)*Aart/Vlng -Qlvr*Alvr/Vlvr +Qspl*Aspr/Vspl;
  dxdt_Abra = Qbra*(Aart/Vlng -Abra/Vbra);
  dxdt_Amsc = Qmsc*(Aart/Vlng -Amsc/Vmsc);
  dxdt_Aspr = Qspl*(Aart/Vlng -Aspr/(Vspl)) +PSspl*(Asps/(Vspl) -Aart/Vlng);
  dxdt_Asps = PSspl*(Aart/Vlng -Asps/(Vspl));
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
  double Cspr = Aspr/(0.25*Vspl);
  double Csps = Aspr/(0.75*Vspl);
  double Cbod = Abod/Vbod;

$CAPTURE
  Cpa Cart Clng Chrt Ckid Clvr Cbra Cmsc Cspr Csps Cbod
  COstd WTstd Qhrt Qkid Qlvr Qbra Qmsc Qspl Qbod Qco
  Vmix Vlng Vhrt Vkid Vlvr Vbra Vmsc Vspl Vbod
'
# Compile code
  brown.mod <- mcode("mouseBROWN", code)
