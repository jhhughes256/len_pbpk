# Define the model parameters and equations
	# Using mrgsolve - analytical solutions
	# Cannot have whitespace before $ blocks
# ------------------------------------------------------------------------------

code <- '
$INIT  // Initial conditions for compartments
	PLA = 0,  // Plasma
	BRA = 0,  // Brain
	LIV = 0,  // Liver
	MSC = 0,  // Muscle
	HRT = 0,  // Heart
	SPL = 0,  // Spleen
	LNG = 0,  // Lung
	KID = 0,  // Kidney
	CAR = 0  // Carcass

$PARAM  // Physiological Parameters (Brown et. al 1997)
  // Regional Blood Flow
	Qco = 1,
	Qbra = 0.033,
	Qlvr = 0.161,
	Qmsc = 0.159,
	Qhrt = 0.066,
	Qspl = 0.01125,  // Imputed from (Davies et. al 1993)
	Qkid = 0.091,
	Qcar = ,  // Remainder of cardiac output

	// Tissue Mass Balance (percentage of mass)
	Vmix = 0.049,
	Vbra = 0.017,
	Vlvr = 0.055,
	Vmsc = 0.384,
	Vhrt = 0.005,
	Vspl = 0.0035,
	Vlng = 0.007,
	Vkid = 0.017,
	Vbod = 0.4625,  // Remainder of volume

	// Physicochemical Parameters

  // Default covariate values for simulation
	WT = 0.02,  // Weight (kg)
	CO = 13.98,  // Cardiac Output (ml/min)

$OMEGA  // Omega covariance block

$OMEGA  // Omega variance

$SIGMA  // Sigma

$MAIN  // Determine covariate values
  // Individual parameter values

$ODE  // Differential equations
	dxdt_DEPOT = -KTR*DEPOT;

$TABLE  // Determine individual predictions
  double IPRE = CENT/V1;
	double DV = IPRE*(1+ERR_PRO);

$CAPTURE  // Capture output
  IPRE DV CL V1 KTR
'
	# Compile the model code
	brown.mod <- mcode("mouseBROWN", code)
