/------------------------------------------------------------------\
|                                                                  |
|             State and Parameter inference on heat shock	   |
|          	model by using particle filter	       	           |
|                                                                  |
|    	      Author: Xin Liu					   |
|             Date: 23/04/2011				           |
|								   |
|                                                                  |
\------------------------------------------------------------------/

This implementation includes five files:
*********************************************************************
There are two main file in this fold, one can simultaneously plot the
results when the programme is running. As a result, this version is 
inefficient as compare to another one without this function.
*********************************************************************


hsMainAllPara:			This is the main file of the implementation which does inference for all parameters unknown.
				It can simultaneously plot the results when the programme is running.

hsMainAllParaWithoutPlot: 	This is the main file without simultaneously plot.

hs_ode: 			The ordinary differential equations for heat shock model.

hs_odeSmc_All_Para: 		The ordinary differential equations for heat shock model by treating parameters as variables.

updatePara_all_para: 		The transition of parameters from which construct the artificial dynamics of parameters.

