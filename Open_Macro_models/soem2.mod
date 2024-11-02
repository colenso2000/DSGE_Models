// This file replicates the model studied in Schmitt-Groh� and Uribe (2003) 
// "Closing small open economy models", JIE, 61, pp. 163-185. 
// This Dynare code solves MODEL 2:
// "A model with debt-elastic interest rate premium"
// Code by Diego Vil�n - Summer 2012
// 

%----------------------------------------------------------------
% 0. Housekeeping:
%----------------------------------------------------------------

close all;
clc;


%----------------------------------------------------------------
% 1. Define variables:
%----------------------------------------------------------------

var  c h y i k a lambda d r riskpremium; 
varexo e;                                     

parameters GAMMA OMEGA RHO SIGMA DELTA PSI_2 ALPHA PHI R_BAR D_BAR BETA;


%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------
GAMMA = 2;			 % Risk aversion parameter 
OMEGA = 1.455;		 % Frisch-elasticity parameter 
RHO = 0.42;		    % Autocorrelation of TFP shock  
SIGMA = 0.0129;	 % Standard deviation of TFP shock 
DELTA = 0.1;		 % Capital depreciation rate 
PSI_2 = 0.000742;  % XXX - from Table #2
ALPHA = 0.32;  	 % Labor share 
PHI = 0.028;  		 % Capital adjustment cost parameter 
R_BAR = 0.04; 		 % World interest rate
D_BAR = 0.7442; 	 % Long-run indebtness level


%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model; 
//1. Evolution of debt stock - eq #4:
d = (1 + exp(r(-1)))*d(-1)- exp(y) + exp(c) + exp(i) + (PHI/2)*(exp(k) - exp(k(-1)))^2; 

//2. Production Function - eq #5:
exp(y) = exp(a)*(exp(k(-1))^ALPHA)*(exp(h)^(1-ALPHA)); 

//3. Law of motion for capital - eq #6
exp(k) = exp(i) + (1-DELTA)*exp(k(-1));  

//4. Euler equation - eq #24:
exp(lambda)= BETA*(1 + exp(r))*exp(lambda(+1));  

//5. Labor FOC - eq #26:  
((exp(c)-((exp(h)^OMEGA)/OMEGA))^(-GAMMA))*(exp(h)^(OMEGA-1))  = exp(lambda)*(1-ALPHA)*exp(y)/exp(h);  

//6. Investment FOC - eq #27:  
exp(lambda)*(1 + PHI*(exp(k) - exp(k(-1)))) = BETA*exp(lambda(+1))*(ALPHA*exp(y(+1))/exp(k) + 1 - DELTA + PHI*(exp(k(+1))-exp(k)));  

//7. Law of motion for TFP - eq #14:
a = RHO*a(-1) + SIGMA * e;  

//8. Country interest rate - eq #23:  
exp(r) = R_BAR + riskpremium; 

// Aid/Bridge equations:
riskpremium = PSI_2*(exp(d - D_BAR)-1); 
exp(lambda) = (exp(c)-((exp(h)^OMEGA)/OMEGA))^(-GAMMA);

end; 


%----------------------------------------------------------------
% 4. Computation:
%----------------------------------------------------------------

// Steady state
steady_state_model; 
   BETA   = 1/(1 + R_BAR); 
   r      = log((1-BETA)/BETA); 
   d      = D_BAR; 
   h      = log(((1-ALPHA)*(ALPHA/(R_BAR + DELTA))^(ALPHA/(1-ALPHA)))^(1/(OMEGA-1))); 
   k      = log(exp(h)/(((R_BAR + DELTA)/ALPHA)^(1/(1-ALPHA)))); 
   y      = log((exp(k)^ALPHA)*(exp(h)^(1 - ALPHA))); 
   i      = log(DELTA*exp(k)); 
   c      = log(exp(y)-exp(i)- R_BAR*d); 
   tb_y   = 1-((exp(c) + exp(i))/exp(y));  
   lambda = log((exp(c)-((exp(h)^OMEGA)/OMEGA))^(-GAMMA)); 
   a      = 0; 
   ca_y   = 0; 
   riskpremium = 0; 
 end; 
 
// Check that this is indeed the steady state
steady;


// Check the Blanchard-Kahn conditions
check;

// Declare a positive temporary technological improvement:
shocks;
var e; stderr 1; 
end;


%----------------------------------------------------------------
% 5. Some Results:
%----------------------------------------------------------------


// Stochastic simulation of the model for 200 periods 
stoch_simul(hp_filter=1600, order=1, irf=200);


