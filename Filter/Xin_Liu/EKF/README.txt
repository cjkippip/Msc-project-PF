/------------------------------------------------------------------\
|                                                                  |
|             State and Parameter inference on heat shock	   |
|          	model by using unscented Kalman filter	       	   |
|                                                                  |
|    	      Author: Xin Liu					   |
|             Date: 23/04/2011				           |
|								   |
|                                                                  |
\------------------------------------------------------------------/

This implementation includes six files.

UKFMainAll: 		The main file for the implementation.
hs_ode:     		The ODEs for heat shock model and it is used for generating synthetic data.
hs_odeUKF_allPara:	The ODEs for heat shock model with associating predicted parameters
hs_F_aug:		Heat shock model by associating augmented states
hs_H_aug:		Observation model by associating augmented states
ukf:			Implementation for UKF algorithm.