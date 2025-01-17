```c
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
  Atubc = 0, // Renal Tubule Cells
  Aurine = 0, // Urine
  Ahrt = 0,  // Heart
  Amsc = 0,  // Muscle
  Abod = 0  // Rest of body

$PARAM
  // Standard Physiological Parameters (Brown et. al 1997)
  WTstd = 25,  // Weight (g)
  COstd = 13.98,  // Cardiac Output (ml/min)

  // Regional Blood Flow (mL/min)
  Qbrastd = 0.4613,
  Qlvrstd = 2.251,
  Qgutstd = 1.5,  // Imputed from Davies et. al
  Qsplstd = 0.1573,  // Imputed from Davies et. al
  PSsplstd = 0.181,
  Qkidstd = 1.272,
  Qhrtstd = 0.9227,
  Qmscstd = 2.223,

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

  // Renal Physiology
  GFRstd = 0.33,  // Glomerular Filtration Rate (ml/min) (Davies et. al)
  PSdiffstd = 0.5,  // Tubular Cell Permeability
  kurinestd = 0.00111,  // Urinary Output (1.6 ml/day)

  // Drug Related Parameters
  fu = 0.6,  // Fraction unbound in Plasma
  fuT = 0.48,  // Fraction unbound in Tissues
  Vmaxt = 6,  // Maximum rate of renal tubular secretion
  kmt = 0.001,  // km of renal tubular secretion

  // Default Covariate Values
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
  double kurine = kurinestd*pow(WT/WTstd,0.75);
  double PSdiff = PSdiffstd*pow(WT/WTstd,0.75);

  // Apparent distribution volumes with allometric scaling
  double Vmix = Vmixstd*pow(WT/WTstd,1);  // No apparent distribution component
  double Vlng = Vlngstd*(fu/fuT)*pow(WT/WTstd,1);
  double Vbra = Vbrastd*(fu/fuT)*pow(WT/WTstd,1);
  double Vlvr = Vlvrstd*(fu/fuT)*pow(WT/WTstd,1);
  double Vgut = Vgutstd*(fu/fuT)*pow(WT/WTstd,1);
  double Vspl = Vsplstd*(fu/fuT)*pow(WT/WTstd,1);
  double Vkid = Vkidstd*(fu/fuT)*pow(WT/WTstd,1);
  double Vhrt = Vhrtstd*(fu/fuT)*pow(WT/WTstd,1);
  double Vmsc = Vmscstd*(fu/fuT)*pow(WT/WTstd,1);
  double Vbod = Vbodstd*(fu/fuT)*pow(WT/WTstd,1);

$ODE
  dxdt_Apa = -Qco*Apa/Vmix +Qbra*Abra/Vbra +Qlvr*Alvr/Vlvr +Qkid*Akid/Vkid +Qhrt*Ahrt/Vhrt +Qmsc*Amsc/Vmsc +Qbod*Abod/Vbod;
  dxdt_Aart = Qco*(Apa/Vmix -Aart/Vlng);
  dxdt_Abra = Qbra*(Aart/Vlng -Abra/Vbra);
  dxdt_Alvr = (Qlvr-(Qspl+Qgut))*Aart/Vlng -Qlvr*Alvr/Vlvr
    +Qspl*Aspr/Vspl +Qgut*Agut/Vgut;
  dxdt_Agut = Qgut*(Aart/Vlng -Agut/Vgut);
  dxdt_Aspr = Qspl*(Aart/Vlng -Aspr/Vspl) +PSspl*(Asps -Aspr);
  dxdt_Asps = PSspl*(Aspr -Asps);
  dxdt_Akid = Qkid*(Aart/Vlng -Akid/Vkid) -GFR*Aart/Vlng +PSdiff*(Atubc -Akid);
  double PStran = Vmaxt/(kmt +Akid/Vkid);
  dxdt_Atubf = GFR*Aart/Vlng +PStran*Atubc -Atubf*kurine;
  dxdt_Atubc = PSdiff*(Akid -Atubc) -PStran*Atubc;
  dxdt_Aurine = Atubf*kurine;
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
  double Ctubf = Atubf;
  double Ctubc = Atubc;
  double Curine = Aurine;
  double Chrt = Ahrt/Vhrt;
  double Cmsc = Amsc/Vmsc;
  double Cbod = Abod/Vbod;

$CAPTURE
  Cpa Cart Cbra Clvr Cgut Cspr Csps Ckid Ctubf Ctubc Curine Chrt Cmsc Cbod
  COstd WTstd Qbra Qlvr Qgut Qspl PSspl Qkid Qhrt Qmsc Qbod Qco
  Vmix Vlng Vbra Vlvr Vspl Vkid Vhrt Vmsc Vbod
```
