// "Solving a Real Business Cycle growth model using Dynare"
// Code by Diego Vilán - Summer 2012
// 

// This Dynare code solves a RBC model with imperfect competition.
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

var y c k i l w r  z;
varexo e;

parameters BETA PSI DELTA ALPHA RHO GAMMA SIGMA EPSILON;


%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------

ALPHA   = 0.33;
BETA    = 0.99;
DELTA   = 0.025;
PSI     = 1.75;
RHO     = 0.95;  
SIGMA   = (0.007/(1-ALPHA));
EPSILON = 10;


%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------
model;
  (1/c) = BETA*(1/c(+1))*(1+r(+1)- DELTA);
  PSI*c/(1-l) = w;
  c + i = y;
  y = (k(-1)^ALPHA)*(exp(z)*l)^(1-ALPHA);
  w = y*((EPSILON-1)/EPSILON)*(1-ALPHA)/l;
  r = y*((EPSILON-1)/EPSILON)*ALPHA/k;
  i = k-(1-DELTA)*k(-1);
  z = RHO*z(-1)+e;
end;


%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------

// Steady state (analytically solved)
initval;
  k = 9;
  c = 0.76;
  l = 0.3;
  w = 2.07;
  r = 0.03;
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
