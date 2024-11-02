// "Solving a Neoclassical growth model using Dynare"
// Code by Diego Vilán - Summer 2012
// 

// This Dynare code solves a neo-classical (determinisitic) growth model.
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

var c k;
varexo z;

parameters ALPHA GAMMA DELTA BETA;


%----------------------------------------------------------------
% 2. Calibration:
%----------------------------------------------------------------

ALPHA = 0.33;
GAMMA = 0.5;
BETA = 0.99;
DELTA = 0.025;


%----------------------------------------------------------------
% 3. The Model:
%----------------------------------------------------------------

model; 
exp(z) * exp(k(-1))^ALPHA = exp(c) + exp(k) - (1 - DELTA)* exp(k(-1)); // Resource constraint
exp(c)^(-GAMMA) = BETA * exp(c(+1))^(-GAMMA)*(ALPHA * exp(z(+1)) * exp(k)^(ALPHA-1) + 1 - DELTA); //Euler equation
end;

// Note: this model is in log-form.

%----------------------------------------------------------------
% 4. Computation:
%----------------------------------------------------------------

// Steady state (analytically solved)
initval;
 k  = log((((1/BETA) + DELTA - 1)/ALPHA)^(1/(ALPHA-1)));
 c  = log(k^ALPHA - DELTA*k);
 z  = log(1); 
end;

// Note: I directly feed the SS into the model, but applying logs


// Check that this is indeed the steady state
steady;


// Check the Blanchard-Kahn conditions
check;


// Declare a positive temporary technological improvement in period 10
shocks;
var z;

// Exact period in which the innovation will take place.
periods 10; 

// periods 10:25; //Aside: You could also define a range instead of 1 period.

values 1.2; // Size of the temprary shock

end;



%----------------------------------------------------------------
% 5. Some Results:
%----------------------------------------------------------------


// Deterministic simulation of the model for 250 periods 

simul(periods = 250);

//Note: for deterministic simulations use simul, rather than stoch_simul

// Display the path of consumption and capital
rplot c;
rplot k;






