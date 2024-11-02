// "Solving a Real Business Cycle growth model using Dynare"
// Code by Diego Vilán - Summer 2012
// 

// This Dynare code solves a basic RBC model featuring log utility in 
// consumption but not in leisure. 
// It features endogenous savings and labor supply.
// The economy is exposed to a 1 period TEMPORARY shock and then simulated.
// See my notes on RBC models for further details.

%----------------------------------------------------------------
% 0. Housekeeping:
%----------------------------------------------------------------

close all;
clc;


%----------------------------------------------------------------
% 1. Define variables:
%----------------------------------------------------------------
var y c k i l z;
varexo e;

parameters BETA PSI DELTA ALPHA RHO GAMMA;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------

ALPHA   = 0.33;
BETA    = 0.99;
DELTA   = 0.025;
PSI     = 1.75;
RHO     = 0.95;  
SIGMA   = (0.007/(1-ALPHA));
GAMMA = 0.5;


%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model; 
 (1/c) = BETA*(1/c(+1))*(1 + ALPHA * y(+1) * k^(-1) - DELTA);
 PSI*c*(l^(1+GAMMA)) = (1-ALPHA)*y;  
  c + i = y;
  y = exp(z)*(k(-1)^ALPHA)*(l^(1-ALPHA));
  i = k -(1-DELTA) * k(-1);
  z = RHO*z(-1) + e;
end;

// 6 equations for 6 unknowns (5 endogenous variables + 1 exogenous one).


%----------------------------------------------------------------
% 4. Computation:
%----------------------------------------------------------------

// Finding the system's steady state:

initval;
 k = 15;
 c = 0.90;
 l = 0.50;
 z = 0; 
 e = 0;
end;

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


%----------------------------------------------------------------
% 6. Some Stats:
%----------------------------------------------------------------

statistic1 = 100*sqrt(diag(oo_.var(1:6,1:6)))./oo_.mean(1:6);
dyntable('Relative standard deviations in %',strvcat('VARIABLE','REL. S.D.'),M_.endo_names(1:6,:),statistic1,10,8,4);
