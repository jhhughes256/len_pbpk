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

$PARAM  // Population parameters (Gerlowski et. al 1983)
  // Regional Blood Flow
	Qco = 1,
	Qbra = 0.033,  // Imputed from (Brown et. al 1997)
	Qlvr = 0.251,
	Qmsc = 0.114,
	Qhrt = 0.0639,
	Qspl = 0.0114,
	Qkid = 0.183,
	Qcar = 0.3437,  // Remainder of cardiac output

	// Tissue Mass Balance (percentage of mass)
	Vmix = 0.045,
	Vbra = 0.017,  // Imputed from (Brown et. al 1997)
	Vlvr = 0.0591,
	Vmsc = 0.455,
	Vhrt = 0.00432,
	Vspl = 0.00455,
	Vlng = 0.00545,
	Vkid = 0.0155,
	Vbod = 0.39408,  // Remainder of volume

	// Covariate effects

	// Default covariate values for simulation
	WT = 0.02,  // Weight (kg)
	CO = 4.38  // Cardiac Output (ml/min)

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
