// "Solving a Neoclassical growth model using Dynare"
// Code by Diego Vilán - Summer 2012
// 

// This Dynare code solves a neo-classical (stochastic) growth model.
// It features endogenous savings decisions, but no endogenous labor supply.
// The economy is exposed to a 1 period TEMPORARY shock and then simulated.
// See my notes on the neo-classical growth model for further details.


%----------------------------------------------------------------
% 0. Housekeeping:
%----------------------------------------------------------------

close all;
clc;


%----------------------------------------------------------------
% 1. Define variables:
%----------------------------------------------------------------

var c k y i z;
varexo e;

parameters ALPHA GAMMA DELTA BETA RHO SIGMA;


%----------------------------------------------------------------
% 2. Calibration:
%----------------------------------------------------------------

ALPHA = 0.33;
GAMMA = 1; //log-utility
BETA = 0.99;
DELTA = 0.025;
SIGMA = 0.01;
RHO = 0.95;


%----------------------------------------------------------------
% 3. The Model:
%----------------------------------------------------------------

model; 
exp(c)^(-GAMMA) = BETA * exp(c(+1))^(-GAMMA)*(ALPHA * exp(z(+1)) * exp(k)^(ALPHA-1) + 1 - DELTA); //Euler equation
exp(y) = exp(z) * exp(k(-1))^(ALPHA);  // Production technology
exp(y) = exp(c) + exp(i); // Resource constraint
exp(k) = exp(i) + (1-DELTA) * exp(k(-1)); // Law of motion of capital
z = RHO * z(-1) + e; // Law of motion of shock
end;

// 5 equations to solve for 5 unknowns (4 endogenous variables, and an exogenous one). 
// Note that in the pdf notes output and investment were substituted out, so there are only 3 equations for 3 unknowns.   
// Note: this model is in log-form.


%----------------------------------------------------------------
% 4. Computation:
%----------------------------------------------------------------

// Steady state (analytically solved)

initval;
  c = log(2.5); 
  k = log(29);
  z = 0;
  i = log(1.5);
  y = log(3);
end;

// Note: The shock does not require logs, because it is originally defined 
// in log terms.


// Check that this is indeed the steady state
steady;


// Check the Blanchard-Kahn conditions
check;


// Declare a positive temporary technological improvement:
shocks;
var e = SIGMA^2;
end;



%----------------------------------------------------------------
% 5. Some Results:
%----------------------------------------------------------------


// Stochastic simulation of the model for 200 periods 

stoch_simul(hp_filter=1600, order=1, irf=200);





