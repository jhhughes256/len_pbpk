# Define the model parameters and equations
	# Using mrgsolve - analytical solutions
	# Cannot have whitespace before $ blocks
# ------------------------------------------------------------------------------

code <- '
$INIT
  // Initial conditions for compartments
	Apa = 0,  // Vascular mixing
	Aart = 0,  // Arterial Blood
	Alng = 0,  // Lung
	Ahrt = 0,  // Heart
	Alvr = 0,  // Liver
	Abod = 0  // Rest of body

$PARAM
  // Standard Physiological Parameters (Brown et. al 1997)
	WTstd = 25,  // Weight (g)
	COstd = 13.98,  // Cardiac Output (ml/min)

  // Regional Blood Flow (mL/min)
	Qhrtstd = 0.9227,
	Qlvrstd = 2.251,

	// Tissue Mass Balance (mL)
	Vmixstd = 1.225,
	Vlngstd = 0.175,
	Vhrtstd = 0.125,
	Vlvrstd = 1.375,

$MAIN
  // Remainder of cardiac output and volume
	double Qbodstd = COstd-(QlvrstdQhrtstd);
	double Vbodstd = WTstd-(Vmixstd+Vlngstd+Vlvrstd+Vhrtstd);

	// Allometric scaling of blood flows, clearances and permeabilities
	double Qhrt = Qhrtstd*pow(WT/WTstd,0.75);
	double Qlvr = Qlvrstd*pow(WT/WTstd,0.75);
	double Qbod = Qbodstd*pow(WT/WTstd,0.75);
	double Qco = Qbra+Qlvr+Qmsc+Qhrt+Qspl+Qkid+Qbod;

	// Apparent distribution volumes with allometric scaling
	double Vmix = Vmixstd*pow(WT/WTstd,1);
	double Vlng = Vlngstd*pow(WT/WTstd,1);
	double Vhrt = Vhrtstd*pow(WT/WTstd,1);
	double Vlvr = Vlvrstd*pow(WT/WTstd,1);
	double Vbod = Vbodstd*pow(WT/WTstd,1);

$ODE
	dxdt_Apa = -Qco*Apa/Vmix +Qhrt*Ahrt/Vhrt +Qkid*Akid/Vkid +Qlvr*Alvr/Vlvr +Qbra*Abra/Vbra +Qmsc*Amsc/Vmsc +Qspl*Aspl/Vspl +Qbod*Abod/Vbod;
	dxdt_Aart = Qco*(Apa/Vmix -Aart/Vlng) +PSlng*(Alng -Aart);
	dxdt_Alng = PSlng*(Aart -Alng);
	dxdt_Ahrt = Qhrt*(Aart/Vlng -Ahrt/Vhrt);
	dxdt_Akid = Qkid*(Aart/Vlng -Akid/Vkid);
	dxdt_Alvr = Qlvr*(Aart/Vlng -Alvr/Vlvr);
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
  Cpa Cart Clng Chrt Ckid Clvr Cbra Cmsc Cspl Cbod COstd WTstd Qhrt
'
	# Compile the model code
	brown.mod <- mcode("mouseBROWN", code)
