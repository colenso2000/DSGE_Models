// This file replicates the model studied in Schmitt-Grohé and Uribe (2003) 
// "Closing small open economy models", JIE, 61, pp. 163-185. 
// This Dynare code solves MODEL 1:
// "A model with an endogenous discount factor"
// Code by Diego Vilán - Summer 2012
// 

%----------------------------------------------------------------
% 0. Housekeeping:
%----------------------------------------------------------------

close all;
clc;


%----------------------------------------------------------------
% 1. Define variables:
%----------------------------------------------------------------
var  c h y i k a d r beta_endo eta lambda util;   
varexo e;                                     

parameters GAMMA OMEGA PSI_1 ALPHA PHI R_BAR DELTA RHO SIGMA D_BAR;


%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------
GAMMA  = 2;     	 % Risk aversion parameter 
OMEGA  = 1.455; 	 % Frisch-elasticity parameter 
PSI_1  = 0.11;       % Elasticity of discount factor - from Table #1 
ALPHA  = 0.32;   	 % Labor share 
PHI    = 0.028;  	 % Capital adjustment cost parameter 
R_BAR    = 0.04; 	 % World interest rate		 
DELTA  = 0.1; 		 % Capital depreciation rate 
RHO    = 0.42;  	 % Autocorrelation of TFP shock  
SIGMA = 0.0129;      % Standard deviation of TFP shock 
D_BAR  = 0.7442;     % Long-run indebtness level


%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------
//1. Evolution of debt stock - eq #4:
d = (1+exp(r(-1)))*d(-1)- exp(y)+exp(c)+exp(i)+(PHI/2)*(exp(k)-exp(k(-1)))^2; 

//2. Production Function - eq #5:
exp(y) = exp(a)*(exp(k(-1))^ALPHA)*(exp(h)^(1-ALPHA)); 

//3. Law of motion for capital - eq #6:
exp(k) = exp(i)+(1-DELTA)*exp(k(-1));  

//4. Euler equation - eq #8:
exp(lambda)= beta_endo*(1+exp(r))*exp(lambda(+1));

//5. Definition marginal utility - eq #9:
exp(lambda)=(exp(c)-((exp(h)^OMEGA)/OMEGA))^(-GAMMA)-exp(eta)*(-PSI_1*(1+exp(c)-OMEGA^(-1)*exp(h)^OMEGA)^(-PSI_1-1));   

//6. LOM for the Lagrange mulitplier on the endogenous discount factor - eq #10:
exp(eta)=-util(+1)+exp(eta(+1))*beta_endo(+1); 

//7. Intratemporal Condition - eq #11:
((exp(c)-(exp(h)^OMEGA)/OMEGA)^(-GAMMA))*(exp(h)^(OMEGA-1)) +  
exp(eta)*(-PSI_1*(1+exp(c)-OMEGA^(-1)*exp(h)^OMEGA)^(-PSI_1-1)*(-exp(h)^(OMEGA-1))) = exp(lambda)*(1-ALPHA)*exp(y)/exp(h);  

//8. Investment Optimality - eq #12:
exp(lambda)*(1+PHI*(exp(k)-exp(k(-1)))) = beta_fun*exp(lambda(+1))*(ALPHA*exp(y(+1))/exp(k)+1-DELTA+PHI*(exp(k(+1))-exp(k)));  

//9. Law of motion for TFP - eq #14:
a = RHO*a(-1)+SIGMA*e;  

//10. Definition endogenous discount factor - p. 168: 
beta_endo =(1+exp(c)-OMEGA^(-1)*exp(h)^OMEGA)^(-PSI_1); 

//11. Country interest rate as World i-rate - eq #13:  
exp(r) = R_BAR; 

// Aid/Bridge equations:
util=(((exp(c)-OMEGA^(-1)*exp(h)^OMEGA)^(1-GAMMA))-1)/(1-GAMMA); 

                                   
end; 
  
%----------------------------------------------------------------
% 4. Computation:
%----------------------------------------------------------------

steady_state_model; 
     r     = log(R_BAR); 
     d     = D_BAR; 
     h     = log(((1-ALPHA)*(ALPHA/(R_BAR+DELTA))^(ALPHA/(1-ALPHA)))^(1/(OMEGA-1))); 
     k     = log(exp(h)/(((R_BAR+DELTA)/ALPHA)^(1/(1-ALPHA)))); 
     y     = log((exp(k)^ALPHA)*(exp(h)^(1-ALPHA))); 
     i     = log(DELTA*exp(k)); 
     c     = log(exp(y)-exp(i)-R_BAR*d);  
     util = (((exp(c)-OMEGA^(-1)*exp(h)^OMEGA)^(1-GAMMA))-1)/(1-GAMMA); 
     PSI_1 = -log(1/(1+R_BAR))/(log((1+exp(c)-OMEGA^(-1)*exp(h)^OMEGA))); 
     beta_fun = (1+exp(c)-OMEGA^(-1)*exp(h)^OMEGA)^(-PSI_1); 
     eta = log(-util/(1-beta_fun)); 
     lambda = log((exp(c)-((exp(h)^OMEGA)/OMEGA))^(-GAMMA)-exp(eta)*(-PSI_1*(1+exp(c)-OMEGA^(-1)*exp(h)^OMEGA)^(-PSI_1-1))); 
     a     = 0; 
 end; 



