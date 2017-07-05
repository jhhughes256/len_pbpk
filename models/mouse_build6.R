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
# Build: Adding kidney model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Model code
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  code <- '
$INIT
  // Initial conditions for compartments
  Apa = 0,  // Vascular mixing
  Aart = 0,  // Arterial Blood
  Abra = 0,  // Brain
  Alvr = 0,  // Liver
  Agut = 0,  // GIT
  Aspr = 0,  // Spleen sinus
  Asps = 0,  // Spleen pulp
  Akid = 0,  // Kidney
  Atubf = 0,  // Renal Tubule Filtrate
  Ahrt = 0,  // Heart
  Amsc = 0,  // Muscle
  Abod = 0  // Rest of body

$PARAM
  // Standard Physiological Parameters (Brown et. al 1997)
  WTstd = 25,  // Weight (g)
  COstd = 13.98,  // Cardiac Output (ml/min)
  GFRstd = 0.33,  // Glomerular Filtration Rate (ml/min) (Davies et. al)

  // Permeability Surface Area (mL/min)
  PSsplstd = 0.0158,

  // Regional Blood Flow (mL/min)
  Qbrastd = 0.4613,
  Qlvrstd = 2.251,
  Qgutstd = 1.5,  // Imputed from Davies et. al
  Qsplstd = 0.1573,  // Imputed from Davies et. al
  Qkidstd = 1.272,
  Qhrtstd = 0.9227,
  Qmscstd = 2.223,
  kurinestd = 0.000694,  // 1 ml/day

  // Tissue Mass Balance (mL)
  Vmixstd = 1.225,
  Vlngstd = 0.175,
  Vbrastd = 0.425,
  Vlvrstd = 1.375,
  Vgutstd = 1.425,
  Vsplstd = 0.0875,
  Vkidstd = 0.425,
  Vhrtstd = 0.125,
  Vmscstd = 9.6,

  // Individual Covariate Values
  WT = 28

$MAIN
  // Remainder of cardiac output and volume
  double Qbodstd = COstd-(Qbrastd+Qlvrstd+Qkidstd+Qhrtstd+Qmscstd);
  double Vbodstd = WTstd-(Vbrastd+Vlvrstd+Vgutstd+Vsplstd+Vkidstd+Vhrtstd+Vmscstd);

  // Allometric scaling of blood flows
  double Qbra = Qbrastd*pow(WT/WTstd,0.75);
  double Qlvr = Qlvrstd*pow(WT/WTstd,0.75);
  double Qgut = Qgutstd*pow(WT/WTstd,0.75);
  double Qspl = Qsplstd*pow(WT/WTstd,0.75);
  double Qkid = Qkidstd*pow(WT/WTstd,0.75);
  double Qhrt = Qhrtstd*pow(WT/WTstd,0.75);
  double Qmsc = Qmscstd*pow(WT/WTstd,0.75);
  double Qbod = Qbodstd*pow(WT/WTstd,0.75);
  double Qco = Qhrt+Qkid+Qlvr+Qbra+Qmsc+Qbod;

  // Allometric scaling for clearances and permeabilities
  double PSspl = PSsplstd*pow(WT/WTstd,0.75);
  double GFR = GFRstd*pow(WT/WTstd,0.75);
  double kurine = kurine*pow(WT/WTstd,0.75);

  // Apparent distribution volumes with allometric scaling
  double Vmix = Vmixstd*pow(WT/WTstd,1);
  double Vlng = Vlngstd*pow(WT/WTstd,1);
  double Vbra = Vbrastd*pow(WT/WTstd,1);
  double Vlvr = Vlvrstd*pow(WT/WTstd,1);
  double Vgut = Vgutstd*pow(WT/WTstd,1);
  double Vspl = Vsplstd*pow(WT/WTstd,1);
  double Vkid = Vkidstd*pow(WT/WTstd,1);
  double Vhrt = Vhrtstd*pow(WT/WTstd,1);
  double Vmsc = Vmscstd*pow(WT/WTstd,1);
  double Vbod = Vbodstd*pow(WT/WTstd,1);

$ODE
  dxdt_Apa = -Qco*Apa/Vmix +Qbra*Abra/Vbra +Qlvr*Alvr/Vlvr +Qkid*Akid/Vkid
    +Qhrt*Ahrt/Vhrt +Qmsc*Amsc/Vmsc +Qbod*Abod/Vbod;
  dxdt_Aart = Qco*(Apa/Vmix -Aart/Vlng);
  dxdt_Abra = Qbra*(Aart/Vlng -Abra/Vbra);
  dxdt_Alvr = (Qlvr-(Qspl+Qgut))*Aart/Vlng -Qlvr*Alvr/Vlvr
    +Qspl*Aspr/Vspl +Qgut*Agut/Vgut;
  dxdt_Agut = Qgut*(Aart/Vlng -Agut/Vgut);
  dxdt_Aspr = Qspl*(Aart/Vlng -Aspr/Vspl) +PSspl*(Asps/Vspl -Aart/Vlng);
  dxdt_Asps = PSspl*(Aart/Vlng -Asps/Vspl);
  dxdt_Akid = Qkid*(Aart/Vlng -Akid/Vkid) -GFR*Aart/Vlng;
  dxdt_Atubf = GFR*Aart/Vlng -Atubf*kurine;
  dxdt_Ahrt = Qhrt*(Aart/Vlng -Ahrt/Vhrt);
  dxdt_Amsc = Qmsc*(Aart/Vlng -Amsc/Vmsc);
  dxdt_Abod = Qbod*(Aart/Vlng -Abod/Vbod);

$TABLE  // Determine individual predictions
  double Cpa = Apa/Vmix;
  double Cart = Aart/Vlng;
  double Cbra = Abra/Vbra;
  double Clvr = Alvr/Vlvr;
  double Cgut = Agut/Vgut;
  double Cspr = Aspr/(0.25*Vspl);
  double Csps = Aspr/(0.75*Vspl);
  double Ckid = Akid/Vkid;
  double Chrt = Ahrt/Vhrt;
  double Cmsc = Amsc/Vmsc;
  double Cbod = Abod/Vbod;

$CAPTURE
  Cpa Cart Cbra Clvr Cgut Cspr Csps Ckid Chrt Cmsc Cbod
  COstd WTstd Qbra Qlvr Qgut Qspl PSspl Qkid Qhrt Qmsc Qbod Qco
  Vmix Vlng Vbra Vlvr Vspl Vkid Vhrt Vmsc Vbod
'
# Compile code
  brown.mod <- mcode("mouseBROWN", code)
