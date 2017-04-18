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

$PARAM  // Population parameters (Davies et. al 1993)
  // Regional Blood Flow
	Qco = 1,
	Qbra = 0.033,  // Imputed from (Brown et. al 1997)
	Qlvr = 0.225,
	Qmsc = 0.11375,
	Qhrt = 0.035,
	Qspl = 0.01125,
	Qkid = 0.1625,
	Qcar = ,  // Remainder of cardiac output

	// Tissue Mass Balance (percentage of mass)
	Vmix = 0.085,
	Vbra = 0.017,  // Imputed from (Brown et. al 1997)
	Vlvr = 0.065,
	Vmsc = 0.5,
	Vhrt = 0.00475,
	Vspl = 0.005,
	Vlng = 0.005,
	Vkid = 0.017,
	Vbod = 0.30125,  // Remainder of volume

	// Covariate effects

	// Default covariate values for simulation
	WT = 0.02,  // Weight (kg)
	CO = 8.0  // Cardiac Output (ml/min)

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
	davies.mod <- mcode("mouseDAVIES", code)
