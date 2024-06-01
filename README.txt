		
About the Direct reduction model
		
	This model is built during the master thesis project "Dynamic modeling and simulation of Direct reduction furnace" in the spring of 2024.
	The model is centered on the 1D mass- and energy balances inside the furnace and the Unreacted shrinking core model for describing the pellet reduction, 
	it also contains tables and expressions for how material constants, equation coefficients etc. changes with temperature and concentration. For full details of the
	model, see the thesis report.

	The folder contains a number of MATLAB scripts and Simulink files that can be used to find steady state, generate datasets, simulate the dynamic behavior
	or regulate the model. All different files will be described below:

	MATLAB Script: calc_lambda_cp_my.m
	
		This script generates datasets and tables that is necessary to the simulation models. The script uses expressions and tabulated values to create datasets 
		that is stored in structs for each element and then saved in the file lambda_cp_my.mat so that the script doesn't need to be run every time.
	
		Example: All necessary data for hydrogen gas can be found th the struct H2. All properties can be found in the different fields of the struct and are tabulated 
		against the temperature, which can be found in H2.T.  

	MATLAB Script: Const.mlx

		This script is used as a initialization script for the SS_solve script. In the file, there are multiple coefficients, constants and reactor-specific parameters that can be 
		modified to change the process. All parameters regarding material properties are added to each material-struct, while the rest of the constants are added to the 
		struct named Var. In some cases, there exists an option to choose between two ways of entering the constant (like choosing between number of pellets or the 
		porosity of the reactor)

		Example: The iron-ore mass flow are a property that can be changed. By changing the ms value in the script, the new value will be stored in the struct. However, the
		mass flow is not used directly in the model and therefore the mass flow is calculated into velocity instead. If the velocity is known, then there exists an option to 
		enter the velocity instead.

	MATLAB Script: coeff.m

		This script calculates the necessary coefficient inside the reactor given values of temperatures and concentrations for the different materials. coeff interpolates
		material properties from the material structs and uses then expressions to find values of the gas mix conductivity, the reaction rate, metallizations grade etc. 

		The script is used together with derivativecalcc in the SS_solve.m script to find the steady state solutions.

	MATLAB Script: derivativecalcc.m

		derivativecalcc calculates the state time-derivatives in the DR-model based on the known PDEs and the upwind discretization of the equations using an input vector X
		containing the temperatures and concentrations in a N-by-5 vector. The full derivation of these matrices can be found in the thesis report. Output gives a value of the 
		time derivative of each inner point of the reactor, stored in a analogous N-by-5 vector.

		The script uses coeff.m to calculate all coefficients and is used in the SS_solve.m script to find steady state solutions

	MATLAB Script: SS_solve.m

		ss

	MATLAB Script: Initialization.mlx

		ss

	Simulink file: Model.slx

		ss

	Simulink file: Model_regulator.slx

		asd