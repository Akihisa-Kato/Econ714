// HW 2 Tangeble and Intangeble Capital
// Firm's decision
// Akihisa Kato
// 

//----------------------------------------------------------------
// 0. Housekeeping
//----------------------------------------------------------------

close all

//----------------------------------------------------------------
// 1. Endogenous variables (18)
//----------------------------------------------------------------
 
var 

// Value (1)
J

// Allocation (8)
kt ki l1 l2 ki1 kt1 s b

// Stochastic Process (2)
z1 z2

// Lagrangean (2)
mmu1 mmu2

// Auxiliary functions (2)
g1 g2

// Envelope condition (3)
Jkt Jki Jb;

//----------------------------------------------------------------
// 2. Exogenous variables (2)
//----------------------------------------------------------------
 
varexo 

ez1 ez2;

//----------------------------------------------------------------
// 3. Parameters
//----------------------------------------------------------------

parameters 

// Technology
llambda1 llambda2 ttheta1 ttheta2 aalpha1 aalpha2 ggamma1 ggamma2
ddeltat ddeltai

// Prices
w ri r 

// Autoregressive stochastic processes
rrhoz1 rrhoz2

// S.D.'s stochastic processes
ssigmaz1 ssigmaz2;


//----------------------------------------------------------------
// 4. Calibration
//----------------------------------------------------------------

// Technology
llambda1    = 1.1;
llambda2    = 0.9;
ttheta1     = 0.5;
ttheta2     = 0.4; 
aalpha1     = 0.3;
aalpha2     = 0.4;
ggamma1     = 0.6;
ggamma2     = 0.5;
ddeltat     = 0.10;
ddeltai     = 0.15;

// Prices
w           = 1;
ri          = 0.025;
r           = 0.03;

// Autoregressive stochastic processes
rrhoz1      = 0.95;
rrhoz2      = 0.90;

// S.D.'s stochastic processes
ssigmaz1    = 0.05;
ssigmaz2    = 0.10;

// Steady state
mmu1_ss = 1;
Jkt_ss = 1+ri;
//Jb_ss = -(1+r);



//----------------------------------------------------------------
// 5. Model 13
//----------------------------------------------------------------

model; 
   
  // 1. Value function
  J = s-w*(l1+l2)+1/(1+ri)*J(+1);  

  // 2. Auxiliary (capital component of tangeble production)
  g1 = ttheta1*kt1^((llambda1-1)/llambda1)+(1-ttheta1)*ki1^((llambda1-1)/llambda1);

  // 3. Auxiliary (capital component of intangeble production)
  g2 = ttheta2*(kt(-1)-kt1)^((llambda2-1)/llambda2)+(1-ttheta2)*(ki(-1)-ki1)^((llambda2-1)/llambda2);
  
  // 4. FOC for kt1
  mmu1*aalpha1*g1^((llambda1*(aalpha1-1)+1)/(llambda1-1))*ttheta1*kt1^(-1/llambda1)*(exp(z1)*l1)^ggamma1 = mmu2*aalpha2*g2^((llambda2*(aalpha2-1)+1)/(llambda2-1))*ttheta2*(kt(-1)-kt1)^(-1/llambda2)*(exp(z2)*l2)^ggamma2;
  
  // 5. FOC for ki1
  mmu1*aalpha1*g1^((llambda1*(aalpha1-1)+1)/(llambda1-1))*(1-ttheta1)*ki1^(-1/llambda1)*(exp(z1)*l1)^ggamma1 = mmu2*aalpha2*g2^((llambda2*(aalpha2-1)+1)/(llambda2-1))*(1-ttheta2)*(ki(-1)-ki1)^(-1/llambda2)*(exp(z2)*l2)^ggamma2;

  // 6. FOC for l1
  w*l1 = mmu1*ggamma1*g1^(aalpha1*llambda1/(llambda1-1))*(exp(z1)*l1)^ggamma1;

  // 7. FOC for l2
  w*l2 = mmu2*ggamma2*g2^(aalpha2*llambda2/(llambda2-1))*(exp(z2)*l2)^ggamma2;

  // 8. FOC for kt'
  Jkt(+1)/(1+ri) = mmu1;
      
  // 9. Envelope theorem for kt
  Jkt = mmu1*(1-ddeltat)+mmu2*aalpha2*g2^((llambda2*(aalpha2-1)+1)/(llambda2-1))*(exp(z2)*l2)^ggamma2*ttheta2*(kt(-1)-kt1)^(-1/llambda2);

  // 10 Production for tangeble
  s+kt-(1-ddeltat)*kt(-1) = g1^(aalpha1*llambda1/(llambda1-1))*(exp(z1)*l1)^ggamma1;

  // 11 Production for intangeble
  ki = g2^(aalpha2*llambda2/(llambda2-1))*(exp(z2)*l2)^ggamma2+(1-ddeltai)*ki(-1);
  
  // 12. Tangeble productivity process 
  z1 = rrhoz1*z1(-1)+ssigmaz1*ez1;

  // 13. Intangeble productivity process 
  z2 = rrhoz2*z2(-1)+ssigmaz2*ez2;

  // 14. FOC ki'
  Jki(+1)/(1+ri) = mmu2;

  // 15. EC for ki
  Jki = mmu2*(1-ddeltai)+mmu2*aalpha2*g2^((llambda2*(aalpha2-1)+1)/(llambda2-1))*(exp(z2)*l2)^ggamma2*(1-ttheta2)*(ki(-1)-ki1)^(-1/llambda2);

  // 16. FOC s
  1 = mmu1;

  // 17. FOC b'
  0.04*(b-b(-1)) = 1+Jb(+1)/(1+ri);
  
  // 18. ET for b
  Jb = -(1+r)+0.04*(b-b(-1))-0.02*(b(-1)-0.2);
  
end;

//----------------------------------------------------------------
// 6. Computation
//----------------------------------------------------------------

initval;
J    		= 0.539877;
kt   		= 0.104443;
ki   		= 0.104617;
l1   		= 0.0443293;
l2   		= 0.00594081;
ki1  		= 0.0833954;
kt1  		= 0.088917;
s    		= 0.0634378;
b           = -0.05;
z1   		= 0;
z2   		= 0;
mmu1 		= 1;
mmu2 		= 0.757152;
g1   		= 0.800183;
g2   		= 1.556;
Jkt  		= 1.025;
Jki  		= 0.776081;
Jb          = -1.025;
end;

shocks;
  var ez1   = 1;
  var ez2   = 1;
end;

steady;


//stoch_simul(hp_filter = 1600, irf = 20, order = 3, pruning, nocorr, nodecomposition, nomoments);
//http://www.dynare.org/dynare-matlab-m2html/matlab/simult_.html#_top
//simult_(oo_.dr.ys,oo_.dr,oo_.exo_simul,3)
//disp_dr(oo_.dr,3,M_.endo_names)
//k_order_pert(oo_.dr,M,options) 
//stochastic_solvers(oo_.dr,0,M_,options_,oo_)

//
// Run Dynare using pruning package
//

// step 1 run 1st order
stoch_simul(irf = 0, order = 1, noprint, nomoments);
f_11 = [oo_.dr.ghx oo_.dr.ghu];
stoch_simul(irf = 0, order = 3, noprint, nomoments, pruning, periods =1000);
// Step 2
optPruning.numSim = 1000;
optPruning.seedNum = 50;
optPruning.orderApp = options_.order;
optPruning.plotIRF = 0;
optPruning.AntitheticShock = 1;
outDynare = RunDynarePruning(optPruning,oo_,M_,f_11);
