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
# Build: Adjustments for IP dosing
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Model code
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  code <- '
$INIT
  // Initial conditions for compartments
  Apa = 0,  // Vascular mixing
  Apo = 0, // Intestinal mixing
  Aart = 0,  // Arterial Blood
  Alvr = 0,  // Liver
  Agut = 0,  // GIT
  Aspr = 0,  // Spleen sinus
  Asps = 0,  // Spleen pulp
  Akid = 0,  // Kidney
  Atubf = 0,  // Renal Tubule Filtrate
  Atubc = 0, // Renal Tubule Cells
  Aurine = 0, // Urine
  Abra = 0,  // Brain
  Ahrt = 0,  // Heart
  Amsc = 0,  // Muscle
  Abod = 0  // Rest of body

$PARAM
  // Standard Physiological Parameters (Davies et. al 1997)
  WTstd = 20,  // Weight (g)
  COstd = 8.0,  // Cardiac Output (ml/min)

  // Regional Blood Flow (mL/min)
  Qlvrstd = 1.8,
  Qgutstd = 1.5,
  Qsplstd = 0.09,
  PSsplstd = 0.009,
  Qkidstd = 1.3,
  Qbrastd = 0.265,  // Imputed from Brown et. al 1997
  Qhrtstd = 0.28,
  Qmscstd = 0.91,

  // Tissue Mass Balance (mL)
  Vmixstd = 1.7,
  Vlngstd = 0.1,
  Vlvrstd = 1.3,
  Vgutstd = 1.5,
  Vsplstd = 0.1,
  Vkidstd = 0.34,
  Vbrastd = 0.36,  // Imputed from Brown et. al 1997
  Vhrtstd = 0.095,
  Vmscstd = 10,

  // Renal Physiology
  CLgfrstd = 0.33,  // Glomerular Filtration Rate (ml/min) (Davies et. al)
  PSdiffstd = 0.5,  // Tubular Cell Permeability
  kurinestd = 0.00111,  // Urinary Output (1.6 ml/day)

  // Drug Related Parameters
  fu = 0.5,  // Fraction unbound in Plasma
  fuT = 0.95,  // Fraction unbound in Tissues
  Vmax = 100,  // Maximum rate of renal tubular secretion
  km = 0.01,  // km of renal tubular secretion
  ka = 0.006,  // Absorption constant

  // Individual Covariate Values
  WT = 28

$MAIN
  // Remainder of cardiac output and volume
  double Qbodstd = COstd-(Qlvrstd+Qkidstd+Qbrastd+Qhrtstd+Qmscstd);
  double Vbodstd = WTstd-(Vlvrstd+Vgutstd+Vsplstd+Vkidstd+Vbrastd+Vhrtstd+Vmscstd);

  // Allometric scaling of blood flows
  double Qlvr = Qlvrstd*pow(WT/WTstd,0.75);
  double Qgut = Qgutstd*pow(WT/WTstd,0.75);
  double Qspl = Qsplstd*pow(WT/WTstd,0.75);
  double Qkid = Qkidstd*pow(WT/WTstd,0.75);
  double Qbra = Qbrastd*pow(WT/WTstd,0.75);
  double Qhrt = Qhrtstd*pow(WT/WTstd,0.75);
  double Qmsc = Qmscstd*pow(WT/WTstd,0.75);
  double Qbod = Qbodstd*pow(WT/WTstd,0.75);
  double Qco = Qbra+Qlvr+Qkid+Qhrt+Qmsc+Qbod;

  // Allometric scaling for clearances and permeabilities
  double PSspl = PSsplstd*pow(WT/WTstd,0.75);
  double CLgfr = CLgfrstd*pow(WT/WTstd,0.75);
  double kurine = kurinestd*pow(WT/WTstd,0.75);
  double PSdiff = PSdiffstd*pow(WT/WTstd,0.75);

  // Apparent distribution volumes with allometric scaling
  double Vmix = Vmixstd*pow(WT/WTstd,1);  // No apparent distribution component
  double Vlng = Vlngstd*(fu/fuT)*pow(WT/WTstd,1);
  double Vlvr = Vlvrstd*(fu/fuT)*pow(WT/WTstd,1);
  double Vgut = Vgutstd*(fu/fuT)*pow(WT/WTstd,1);
  double Vspl = Vsplstd*(fu/fuT)*pow(WT/WTstd,1);
  double Vspr = 0.25*Vspl;
  double Vsps = 0.75*Vspl;
  double Vkid = Vkidstd*(fu/fuT)*pow(WT/WTstd,1);
  double Vbra = Vbrastd*(fu/fuT)*pow(WT/WTstd,1);
  double Vhrt = Vhrtstd*(fu/fuT)*pow(WT/WTstd,1);
  double Vmsc = Vmscstd*(fu/fuT)*pow(WT/WTstd,1);
  double Vbod = Vbodstd*(fu/fuT)*pow(WT/WTstd,1);

$ODE
  dxdt_Apa = -Qco*Apa/Vmix +Qbra*Abra/Vbra +Qlvr*Alvr/Vlvr +Qkid*Akid/Vkid +Qhrt*Ahrt/Vhrt +Qmsc*Amsc/Vmsc +Qbod*Abod/Vbod;
  dxdt_Apo = -ka*Apo;
  dxdt_Aart = Qco*(Apa/Vmix -Aart/Vlng);
  dxdt_Alvr = (Qlvr-(Qspl+Qgut))*Aart/Vlng -Qlvr*Alvr/Vlvr +Qspl*Aspr/Vspl +Qgut*Agut/Vgut;
  dxdt_Agut = Qgut*(Aart/Vlng -Agut/Vgut) +ka*Apo;
  dxdt_Aspr = Qspl*(Aart/Vlng -Aspr/Vspl) -PSspl*fu*(Aspr/Vspr -Asps/Vsps);
  dxdt_Asps = PSspl*fuT*(Aspr/Vspr -Asps/Vsps);
  dxdt_Akid = Qkid*(Aart/Vlng -Akid/Vkid) -CLgfr*fu*Aart/Vlng -PSdiff*fu*(Akid -Atubc);
  double ktran = Vmax/(km +Akid/Vkid);
  dxdt_Atubf = CLgfr*fu*Aart/Vlng +ktran*Atubc -Atubf*kurine;
  dxdt_Atubc = PSdiff*fu*(Akid -Atubc) -ktran*Atubc;
  dxdt_Aurine = Atubf*kurine;
  dxdt_Abra = Qbra*(Aart/Vlng -Abra/Vbra);
  dxdt_Ahrt = Qhrt*(Aart/Vlng -Ahrt/Vhrt);
  dxdt_Amsc = Qmsc*(Aart/Vlng -Amsc/Vmsc);
  dxdt_Abod = Qbod*(Aart/Vlng -Abod/Vbod);

$TABLE  // Determine individual predictions
  double Cpa = Apa/Vmix;
  double Cpo = Apo;
  double Cart = Aart/Vlng;
  double Clvr = Alvr/Vlvr;
  double Cgut = Agut/Vgut;
  double Cspr = Aspr/Vspr;
  double Csps = Asps/Vsps;
  double Ckid = Akid/Vkid;
  double Ctubf = Atubf;
  double Ctubc = Atubc;
  double Curine = Aurine;
  double Cbra = Abra/Vbra;
  double Chrt = Ahrt/Vhrt;
  double Cmsc = Amsc/Vmsc;
  double Cbod = Abod/Vbod;

$CAPTURE
  Cpa Cpo Cart Clvr Cgut Cspr Csps Ckid Ctubf Ctubc Curine Cbra Chrt Cmsc Cbod
  COstd WTstd Qlvr Qgut Qspl PSspl Qkid Qbra Qhrt Qmsc Qbod Qco
  Vmix Vlng Vlvr Vspl Vkid Vbra Vhrt Vmsc Vbod
'
# Compile code
  davies.mod <- mcode("mouseDAVIES", code)
