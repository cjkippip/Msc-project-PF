/------------------------------------------------------------------\
|                                                                  |
|             State and Parameter inference on heat shock	   |
|          	model by using extended Kalman filter	       	   |
|                                                                  |
|    	      Author: Xin Liu					   |
|             Date: 23/04/2011				           |
|								   |
|                                                                  |
\------------------------------------------------------------------/

This implementation includes four files.

hsEKFMain_AllPara: 		The main file for the implementation.
hs_odeEKF_all_Para: 		The ODEs for heat shock model and it is used for generating synthetic data.
Jacobian_HeatShock_allparams: 	Jacobian matrix for heat shock model.
dreof: 				The differential Lyapunov equation, using for covariance matrix transition.
DifferentialLyapunov: 		Integration of differential Lyapunov equation and tranforms the size of solution 
		      		as Jacobiam matrix of process model.