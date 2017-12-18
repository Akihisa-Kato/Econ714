% Code for Econ714 HW2 Q1
% Firm's decision


%% Preliminary
clear all;
close all; clc;
tic

fsize = 5;
fsize1 = 3;
lwidth = 1.5;
path_Main=pwd;
cd ../.
path_Figure=[pwd '/Main/Figure'];

addpath(path_Main);
addpath(genpath(path_Figure));

cd(path_Main);


%% 1. Prameters
global llambda1 llambda2 ttheta1 ttheta2 aalpha1 aalpha2 ggamma1 ggamma2 ddeltat ddeltai w ri r nVariable ...
    vProductivityT vProductivityI mTransitionT mTransitionI nMaxOrder nGridProductivityI nGridProductivityT ...
    nIndexJ nIndexkt nIndexl1 nIndexl2 nIndexki2 nIndexkt2 nIndexb nIndexki...
    nTangebleCapitalMax nTangebleCapitalMin nIntangebleCapitalMax nIntangebleCapitalMin nDebtMax nDebtMin mSmolyakGrid...
    nGridSmolyak vChebychev_kt vChebychev_ki vChebychev_b nGridProductivity

% Technology
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

% Prices
w           = 1;
ri          = 0.025;
r           = 0.03;


% Approximation Parameters
nState = 3; % # of continuous state variables
nApprox = 2; % level of apprximation (5 = q = nState + nApprox)


vProductivityT = [log(1.05), log(0.95)];
vProductivityI = [log(1.10), log(0.90)];
nGridProductivityI = length(vProductivityI);
nGridProductivityT = length(vProductivityT);


mTransitionT = [0.95 0.05;
                0.05 0.95];
            
mTransitionI = [0.90 0.10;
                0.10 0.90];
            
nGridProductivity = nGridProductivityI*nGridProductivityT;


%% 2. Steady State 

% Based on Dynare's result from Q2
J_ss    		= 0.539877;
kt_ss   		= 0.104443;
ki_ss   		= 0.104617;
l1_ss   		= 0.0443293;
l2_ss   		= 0.00594081;
ki1_ss  		= 0.0833954;
kt1_ss  		= 0.088917;
ki2_ss  		= ki_ss-ki1_ss;
kt2_ss  		= kt_ss-kt1_ss;

s_ss    		= 0.0634378;
b_ss           = -0.05;
z1_ss   		= 0;
z2_ss   		= 0;
mmu1_ss 		= 1;
mmu2_ss 		= 0.757152;
g1_ss   		= 0.800183; % Inside of CES of tangeble capital production
g2_ss   		= 1.556;    % Inside of CES of intangeble capital production
Jkt_ss  		= 1.025;
Jki_ss  		= 0.776081;
Jb_ss          = -1.025;

%vSteadyState = [J_ss kt_ss ki_ss l1_ss l2_ss ki1_ss kt1_ss s_ss b_ss mmu2_ss Jkt_ss Jki_ss Jb_ss]';
vSteadyState = [J_ss kt_ss kt2_ss ki2_ss l1_ss l2_ss b_ss ki_ss]';

nIndexJ = 1; nIndexkt = 2; nIndexkt2 = 3; nIndexki2 = 4; nIndexl1 = 5; nIndexl2 = 6;
nIndexb = 7; nIndexki = 8;


nVariable = length(vSteadyState);
nTangebleCapitalMax = 1.3*kt_ss;
nTangebleCapitalMin = 0.7*kt_ss;

nIntangebleCapitalMax = 1.3*ki_ss;
nIntangebleCapitalMin = 0.7*ki_ss;

nDebtMin = 1.3*b_ss;
nDebtMax = 0.7*b_ss;

%% 3. Sparse Grids

mIndex = Smolyak_Elem_Isotrop(nState,nApprox);
mSmolyakGrid = Smolyak_Grid(nState,nApprox,mIndex);
nGridSmolyak = size(mSmolyakGrid,1);
%% 4. Smolyak algorithm + Chebychev polinomials
% Note that G(5,3) implies we will find 25 unknown coefficients for one decision rules or functions
% and we will approximate 8 equatins (J, kt', kt1, ki1, l1, l2, b', mmu2)
% conditions. This number includes Lagrangian multipliers and Envelope
% conditions. Moreover, we have 4 shock states.
% As a result, the number of unknown coefficients is 25*13*4.




% I will define mTheta as a three dimensional object, mTheta = (# of
% conditinos)*G(5,3)*(# of states)


% 1st dimension is a variables in the model.
% [J kt kt2 ki2 l1 l2 b mmu2]

% 2nd dimension is the coefficient of theta on a particular
% function. coefficient will be ordered as a lexoicographic order.

% 3rd dimension is a state
% [both high, z1 high, z2 high, both low]


% preliminary
nMaxOrder = 2^(3-1)+1; % we use up to m_3


% Initial Guess for Chebyshev coefficients
mTheta_guess = zeros(nVariable,nGridSmolyak,nGridProductivityT*nGridProductivityI);
mTheta_guess(:,1,1) = vSteadyState;
mTheta_guess(:,1,2) = vSteadyState;
mTheta_guess(:,1,3) = vSteadyState;
mTheta_guess(:,1,4) = vSteadyState;

vChebychev_kt = ones(nGridSmolyak,nMaxOrder);
vChebychev_ki = ones(nGridSmolyak,nMaxOrder);
vChebychev_b  = ones(nGridSmolyak,nMaxOrder);
vChebychev_kt(:,2) = mSmolyakGrid(:,1);
vChebychev_ki(:,2) = mSmolyakGrid(:,2);
vChebychev_b(:,2)  = mSmolyakGrid(:,3);

for i1=3:nMaxOrder
    vChebychev_kt(:,i1) = 2*vChebychev_kt(:,2).*vChebychev_kt(:,i1-1)-vChebychev_kt(:,i1-2);
    vChebychev_ki(:,i1) = 2*vChebychev_ki(:,2).*vChebychev_ki(:,i1-1)-vChebychev_ki(:,i1-2);
    vChebychev_b(:,i1)  = 2*vChebychev_b(:,2).*vChebychev_b(:,i1-1)-vChebychev_b(:,i1-2);
end





%% Solve for Chebyshev coefficients
%options=optimset('Display','Iter','TolFun',10^(-15),'TolX',10^(-15));
mTheta = Residual_Function(mTheta_guess);
Q1.Coeff = mTheta;
%save('Result.mat','Q1','-append');
toc
%% Plot Functions


% TangebleCapiital
nGridTangebleCapital = 150;
vGridTangebleCapital = linspace(nTangebleCapitalMin,nTangebleCapitalMax,nGridTangebleCapital);

% scale down into [-1,1] region

vGridTangebleCapitalScaledown = 2*(vGridTangebleCapital-nTangebleCapitalMin)/(nTangebleCapitalMax-nTangebleCapitalMin)-1;


for nTangebleCapital = 1:nGridTangebleCapital
    % compute Chebychev pol
    vChebychev_kt_plot = ones(1,nMaxOrder);
    vChebychev_ki_plot = ones(1,nMaxOrder);
    vChebychev_b_plot  = ones(1,nMaxOrder);
    vChebychev_kt_plot(2) = vGridTangebleCapitalScaledown(nTangebleCapital);
    vChebychev_ki_plot(2) = 0;
    vChebychev_b_plot(2)  = 0;
    for i1=3:nMaxOrder
        vChebychev_kt_plot(i1) = 2*vChebychev_kt_plot(2)*vChebychev_kt_plot(i1-1)-vChebychev_kt_plot(i1-2);
        vChebychev_ki_plot(i1) = 2*vChebychev_ki_plot(2)*vChebychev_ki_plot(i1-1)-vChebychev_ki_plot(i1-2);
        vChebychev_b_plot(i1)  = 2*vChebychev_b_plot(2)*vChebychev_b_plot(i1-1)-vChebychev_b_plot(i1-2);
    end
    
    vChebychevToday = [1; vChebychev_kt_plot(2:nMaxOrder)'; vChebychev_ki_plot(2:nMaxOrder)'; vChebychev_b_plot(2:nMaxOrder)'; % const and no cross terms (13)
                       vChebychev_kt_plot(2)*vChebychev_ki_plot(2); vChebychev_kt_plot(2)*vChebychev_b_plot(2); vChebychev_ki_plot(2)*vChebychev_b_plot(2); % cross 1*1 (3)
                       vChebychev_kt_plot(2)*vChebychev_ki_plot(3); vChebychev_kt_plot(2)*vChebychev_b_plot(3); vChebychev_ki_plot(2)*vChebychev_kt_plot(3); 
                       vChebychev_ki_plot(2)*vChebychev_b_plot(3);  vChebychev_b_plot(2)*vChebychev_kt_plot(3); vChebychev_b_plot(2)*vChebychev_ki_plot(3); % cross 1*2 (6)
                       vChebychev_kt_plot(3)*vChebychev_ki_plot(3); vChebychev_kt_plot(3)*vChebychev_b_plot(3); vChebychev_ki_plot(3)*vChebychev_b_plot(3)];% cross 2*2 (3)
    for nProdutivity = 1:nGridProductivity
       Q1.functions.vValueFunction.TangebleCapital(nTangebleCapital,nProdutivity) = vec(mTheta(nIndexJ,:,nProdutivity))'*vChebychevToday;
       Q1.functions.vTangeblePolicyFunction.TangebleCapital(nTangebleCapital,nProdutivity) = vec(mTheta(nIndexkt,:,nProdutivity))'*vChebychevToday;
       Q1.functions.vIntangeblePolicyFunction.TangebleCapital(nTangebleCapital,nProdutivity) = vec(mTheta(nIndexki,:,nProdutivity))'*vChebychevToday;
       Q1.functions.vDebtPolicyFunction.TangebleCapital(nTangebleCapital,nProdutivity) = vec(mTheta(nIndexb,:,nProdutivity))'*vChebychevToday;
       Q1.functions.vLaborOnePolicyFunction.TangebleCapital(nTangebleCapital,nProdutivity) = vec(mTheta(nIndexl1,:,nProdutivity))'*vChebychevToday;
       Q1.functions.vLaborTwoPolicyFunction.TangebleCapital(nTangebleCapital,nProdutivity) = vec(mTheta(nIndexl2,:,nProdutivity))'*vChebychevToday;
       Q1.functions.vTangebleCapitalTwoPolicyFunction.TangebleCapital(nTangebleCapital,nProdutivity) = vec(mTheta(nIndexkt2,:,nProdutivity))'*vChebychevToday;
       Q1.functions.vIntangebleCapitalTwoPolicyFunction.TangebleCapital(nTangebleCapital,nProdutivity) = vec(mTheta(nIndexki2,:,nProdutivity))'*vChebychevToday;
       
    end
end

Q1.functions.vTangebleCapitalOnePolicyFunction.TangebleCapital = vGridTangebleCapital'-Q1.functions.vTangebleCapitalTwoPolicyFunction.TangebleCapital;
Q1.functions.vIntangebleCapitalOnePolicyFunction.TangebleCapital = ki_ss-Q1.functions.vIntangebleCapitalTwoPolicyFunction.TangebleCapital;


% IntangebleCapiital
nGridIntangebleCapital = 150;
vGridIntangebleCapital = linspace(nIntangebleCapitalMin,nIntangebleCapitalMax,nGridIntangebleCapital);

% scale down into [-1,1] region

vGridIntangebleCapitalScaledown = 2*(vGridIntangebleCapital-nIntangebleCapitalMin)/(nIntangebleCapitalMax-nIntangebleCapitalMin)-1;


for nIntangebleCapital = 1:nGridIntangebleCapital
    % compute Chebychev pol
    vChebychev_kt_plot = ones(1,nMaxOrder);
    vChebychev_ki_plot = ones(1,nMaxOrder);
    vChebychev_b_plot  = ones(1,nMaxOrder);
    vChebychev_kt_plot(2) = 0;
    vChebychev_ki_plot(2) = vGridIntangebleCapitalScaledown(nIntangebleCapital);
    vChebychev_b_plot(2)  = 0;
    for i1=3:nMaxOrder
        vChebychev_kt_plot(i1) = 2*vChebychev_kt_plot(2)*vChebychev_kt_plot(i1-1)-vChebychev_kt_plot(i1-2);
        vChebychev_ki_plot(i1) = 2*vChebychev_ki_plot(2)*vChebychev_ki_plot(i1-1)-vChebychev_ki_plot(i1-2);
        vChebychev_b_plot(i1)  = 2*vChebychev_b_plot(2)*vChebychev_b_plot(i1-1)-vChebychev_b_plot(i1-2);
    end
    
    vChebychevToday = [1; vChebychev_kt_plot(2:nMaxOrder)'; vChebychev_ki_plot(2:nMaxOrder)'; vChebychev_b_plot(2:nMaxOrder)'; % const and no cross terms (13)
                       vChebychev_kt_plot(2)*vChebychev_ki_plot(2); vChebychev_kt_plot(2)*vChebychev_b_plot(2); vChebychev_ki_plot(2)*vChebychev_b_plot(2); % cross 1*1 (3)
                       vChebychev_kt_plot(2)*vChebychev_ki_plot(3); vChebychev_kt_plot(2)*vChebychev_b_plot(3); vChebychev_ki_plot(2)*vChebychev_kt_plot(3); 
                       vChebychev_ki_plot(2)*vChebychev_b_plot(3);  vChebychev_b_plot(2)*vChebychev_kt_plot(3); vChebychev_b_plot(2)*vChebychev_ki_plot(3); % cross 1*2 (6)
                       vChebychev_kt_plot(3)*vChebychev_ki_plot(3); vChebychev_kt_plot(3)*vChebychev_b_plot(3); vChebychev_ki_plot(3)*vChebychev_b_plot(3)];% cross 2*2 (3)
    for nProdutivity = 1:nGridProductivity
       Q1.functions.vValueFunction.IntangebleCapital(nIntangebleCapital,nProdutivity) = vec(mTheta(nIndexJ,:,nProdutivity))'*vChebychevToday;
       Q1.functions.vTangeblePolicyFunction.IntangebleCapital(nIntangebleCapital,nProdutivity) = vec(mTheta(nIndexkt,:,nProdutivity))'*vChebychevToday;
       Q1.functions.vIntangeblePolicyFunction.IntangebleCapital(nIntangebleCapital,nProdutivity) = vec(mTheta(nIndexki,:,nProdutivity))'*vChebychevToday;
       Q1.functions.vDebtPolicyFunction.IntangebleCapital(nIntangebleCapital,nProdutivity) = vec(mTheta(nIndexb,:,nProdutivity))'*vChebychevToday;
       Q1.functions.vLaborOnePolicyFunction.IntangebleCapital(nIntangebleCapital,nProdutivity) = vec(mTheta(nIndexl1,:,nProdutivity))'*vChebychevToday;
       Q1.functions.vLaborTwoPolicyFunction.IntangebleCapital(nIntangebleCapital,nProdutivity) = vec(mTheta(nIndexl2,:,nProdutivity))'*vChebychevToday;
       Q1.functions.vTangebleCapitalTwoPolicyFunction.IntangebleCapital(nIntangebleCapital,nProdutivity) = vec(mTheta(nIndexkt2,:,nProdutivity))'*vChebychevToday;
       Q1.functions.vIntangebleCapitalTwoPolicyFunction.IntangebleCapital(nIntangebleCapital,nProdutivity) = vec(mTheta(nIndexki2,:,nProdutivity))'*vChebychevToday;
       
    end
end

Q1.functions.vTangebleCapitalOnePolicyFunction.IntangebleCapital = kt_ss-Q1.functions.vTangebleCapitalTwoPolicyFunction.IntangebleCapital;
Q1.functions.vIntangebleCapitalOnePolicyFunction.IntangebleCapital = vGridIntangebleCapital'-Q1.functions.vIntangebleCapitalTwoPolicyFunction.IntangebleCapital;


% Debt
nGridDebt = 150;
vGridDebt = linspace(nDebtMin,nDebtMax,nGridDebt);

% scale down into [-1,1] region

vGridDebtScaledown = 2*(vGridDebt-nDebtMin)/(nDebtMax-nDebtMin)-1;


for nDebt = 1:nGridDebt
    % compute Chebychev pol
    vChebychev_kt_plot = ones(1,nMaxOrder);
    vChebychev_ki_plot = ones(1,nMaxOrder);
    vChebychev_b_plot  = ones(1,nMaxOrder);
    vChebychev_kt_plot(2) = 0;
    vChebychev_ki_plot(2) = 0;
    vChebychev_b_plot(2)  = vGridDebtScaledown(nDebt);
    for i1=3:nMaxOrder
        vChebychev_kt_plot(i1) = 2*vChebychev_kt_plot(2)*vChebychev_kt_plot(i1-1)-vChebychev_kt_plot(i1-2);
        vChebychev_ki_plot(i1) = 2*vChebychev_ki_plot(2)*vChebychev_ki_plot(i1-1)-vChebychev_ki_plot(i1-2);
        vChebychev_b_plot(i1)  = 2*vChebychev_b_plot(2)*vChebychev_b_plot(i1-1)-vChebychev_b_plot(i1-2);
    end
    
    vChebychevToday = [1; vChebychev_kt_plot(2:nMaxOrder)'; vChebychev_ki_plot(2:nMaxOrder)'; vChebychev_b_plot(2:nMaxOrder)'; % const and no cross terms (13)
                       vChebychev_kt_plot(2)*vChebychev_ki_plot(2); vChebychev_kt_plot(2)*vChebychev_b_plot(2); vChebychev_ki_plot(2)*vChebychev_b_plot(2); % cross 1*1 (3)
                       vChebychev_kt_plot(2)*vChebychev_ki_plot(3); vChebychev_kt_plot(2)*vChebychev_b_plot(3); vChebychev_ki_plot(2)*vChebychev_kt_plot(3); 
                       vChebychev_ki_plot(2)*vChebychev_b_plot(3);  vChebychev_b_plot(2)*vChebychev_kt_plot(3); vChebychev_b_plot(2)*vChebychev_ki_plot(3); % cross 1*2 (6)
                       vChebychev_kt_plot(3)*vChebychev_ki_plot(3); vChebychev_kt_plot(3)*vChebychev_b_plot(3); vChebychev_ki_plot(3)*vChebychev_b_plot(3)];% cross 2*2 (3)
    for nProdutivity = 1:nGridProductivity
       Q1.functions.vValueFunction.Debt(nDebt,nProdutivity) = vec(mTheta(nIndexJ,:,nProdutivity))'*vChebychevToday;
       Q1.functions.vTangeblePolicyFunction.Debt(nDebt,nProdutivity) = vec(mTheta(nIndexkt,:,nProdutivity))'*vChebychevToday;
       Q1.functions.vIntangeblePolicyFunction.Debt(nDebt,nProdutivity) = vec(mTheta(nIndexki,:,nProdutivity))'*vChebychevToday;
       Q1.functions.vDebtPolicyFunction.Debt(nDebt,nProdutivity) = vec(mTheta(nIndexb,:,nProdutivity))'*vChebychevToday;
       Q1.functions.vLaborOnePolicyFunction.Debt(nDebt,nProdutivity) = vec(mTheta(nIndexl1,:,nProdutivity))'*vChebychevToday;
       Q1.functions.vLaborTwoPolicyFunction.Debt(nDebt,nProdutivity) = vec(mTheta(nIndexl2,:,nProdutivity))'*vChebychevToday;
       Q1.functions.vTangebleCapitalTwoPolicyFunction.Debt(nDebt,nProdutivity) = vec(mTheta(nIndexkt2,:,nProdutivity))'*vChebychevToday;
       Q1.functions.vIntangebleCapitalTwoPolicyFunction.Debt(nDebt,nProdutivity) = vec(mTheta(nIndexki2,:,nProdutivity))'*vChebychevToday;
       
    end
end

Q1.functions.vTangebleCapitalOnePolicyFunction.Debt = kt_ss-Q1.functions.vTangebleCapitalTwoPolicyFunction.Debt;
Q1.functions.vIntangebleCapitalOnePolicyFunction.Debt = ki_ss-Q1.functions.vIntangebleCapitalTwoPolicyFunction.Debt;


save('Result.mat','Q1','-append');

%%
figure1 = figure(1);
name1 = strcat('ValueFunction_Chebychev.pdf');

subplot(2,2,1)
plot(vGridTangebleCapital',Q1.functions.vValueFunction.TangebleCapital','LineWidth',lwidth)
axis tight;
legend('both high','z1 high','z2 high','both low')
set(gca,'FontSize',fsize)
title('Value Function')
xlabel('Tangeble Capital','FontSize',fsize)
set(legend,'FontSize',fsize1)
set(legend,'Position',[0.22,0.85,0.01,0.01])
set(legend,'Box','off')


subplot(2,2,2)
plot(vGridIntangebleCapital',Q1.functions.vValueFunction.IntangebleCapital','LineWidth',lwidth)
axis tight;
set(gca,'FontSize',fsize)
title('Value Function')
xlabel('Intangeble Capital','FontSize',fsize)

subplot(2,2,3)
plot(vGridDebt',Q1.functions.vValueFunction.Debt','LineWidth',lwidth)
axis tight;
set(gca,'FontSize',fsize)
title('Value Function')
xlabel('Debt','FontSize',fsize)

cd(path_Figure);
set(gcf,'PaperOrientation','landscape','PaperPosition',[0 0 16 12])
set(gcf, 'PaperSize', [16 12]);
print(figure1,'-dpdf',name1)
cd(path_Main);

figure2 = figure(2);
name2 = strcat('TangeblePolicyFunction_Chebychev.pdf');

subplot(2,2,1)
plot(vGridTangebleCapital',Q1.functions.vTangeblePolicyFunction.TangebleCapital','LineWidth',lwidth)
axis tight;
legend('both high','z1 high','z2 high','both low')
set(gca,'FontSize',fsize)
title('Tangeble Policy Function')
xlabel('Tangeble Capital','FontSize',fsize)
set(legend,'FontSize',fsize1)
set(legend,'Position',[0.22,0.85,0.01,0.01])
set(legend,'Box','off')


subplot(2,2,2)
plot(vGridIntangebleCapital',Q1.functions.vTangeblePolicyFunction.IntangebleCapital','LineWidth',lwidth)
axis tight;
set(gca,'FontSize',fsize)
title('Tangeble Policy Function')
xlabel('Intangeble Capital','FontSize',fsize)

subplot(2,2,3)
plot(vGridDebt',Q1.functions.vTangeblePolicyFunction.Debt','LineWidth',lwidth)
axis tight;
set(gca,'FontSize',fsize)
title('Tangeble Policy Function')
xlabel('Intangeble Capital','FontSize',fsize)

cd(path_Figure);
set(gcf,'PaperOrientation','landscape','PaperPosition',[0 0 16 12])
set(gcf, 'PaperSize', [16 12]);
print(figure2,'-dpdf',name2)
cd(path_Main);

figure3 = figure(3);
name3 = strcat('IntangeblePolicyFunction_Chebychev.pdf');

subplot(2,2,1)
plot(vGridTangebleCapital',Q1.functions.vIntangeblePolicyFunction.TangebleCapital','LineWidth',lwidth)
axis tight;
legend('both high','z1 high','z2 high','both low')
set(gca,'FontSize',fsize)
title('Intangeble Pol Fun')
xlabel('Tangeble Capital','FontSize',fsize)
set(legend,'FontSize',fsize1)
set(legend,'Position',[0.22,0.85,0.01,0.01])
set(legend,'Box','off')


subplot(2,2,2)
plot(vGridIntangebleCapital',Q1.functions.vIntangeblePolicyFunction.IntangebleCapital','LineWidth',lwidth)
axis tight;
set(gca,'FontSize',fsize)
title('Intangeble Pol Fun')
xlabel('Intangeble Capital','FontSize',fsize)

subplot(2,2,3)
plot(vGridDebt',Q1.functions.vIntangeblePolicyFunction.Debt','LineWidth',lwidth)
axis tight;
set(gca,'FontSize',fsize)
title('Intangeble Pol Fun')
xlabel('Debt','FontSize',fsize)

cd(path_Figure);
set(gcf,'PaperOrientation','landscape','PaperPosition',[0 0 16 12])
set(gcf, 'PaperSize', [16 12]);
print(figure3,'-dpdf',name3)
cd(path_Main);

figure4 = figure(4);
name4 = strcat('LaborOnePolicyFunction_Chebychev.pdf');

subplot(2,2,1)
plot(vGridTangebleCapital',Q1.functions.vLaborOnePolicyFunction.TangebleCapital','LineWidth',lwidth)
axis tight;
legend('both high','z1 high','z2 high','both low')
set(gca,'FontSize',fsize)
title('Labor 1 Pol Fun')
xlabel('Tangeble Capital','FontSize',fsize)
set(legend,'FontSize',fsize1)
set(legend,'Position',[0.22,0.85,0.01,0.01])
set(legend,'Box','off')


subplot(2,2,2)
plot(vGridIntangebleCapital',Q1.functions.vLaborOnePolicyFunction.IntangebleCapital','LineWidth',lwidth)
axis tight;
set(gca,'FontSize',fsize)
title('Labor 1 Pol Fun')
xlabel('Intangeble Capital','FontSize',fsize)

subplot(2,2,3)
plot(vGridDebt',Q1.functions.vLaborOnePolicyFunction.Debt','LineWidth',lwidth)
axis tight;
set(gca,'FontSize',fsize)
title('Labor 1 Pol Fun')
xlabel('Debt','FontSize',fsize)

cd(path_Figure);
set(gcf,'PaperOrientation','landscape','PaperPosition',[0 0 16 12])
set(gcf, 'PaperSize', [16 12]);
print(figure4,'-dpdf',name4)
cd(path_Main);


figure5 = figure(5);
name5 = strcat('LaborTwoPolicyFunction_Chebychev.pdf');

subplot(2,2,1)
plot(vGridTangebleCapital',Q1.functions.vLaborTwoPolicyFunction.TangebleCapital','LineWidth',lwidth)
axis tight;
legend('both high','z1 high','z2 high','both low')
set(gca,'FontSize',fsize)
title('Labor 2 Pol Fun')
xlabel('Tangeble Capital','FontSize',fsize)
set(legend,'FontSize',fsize1)
set(legend,'Position',[0.22,0.85,0.01,0.01])
set(legend,'Box','off')


subplot(2,2,2)
plot(vGridIntangebleCapital',Q1.functions.vLaborTwoPolicyFunction.IntangebleCapital','LineWidth',lwidth)
axis tight;
set(gca,'FontSize',fsize)
title('Labor 2 Pol Fun')
xlabel('Intangeble Capital','FontSize',fsize)

subplot(2,2,3)
plot(vGridDebt',Q1.functions.vLaborTwoPolicyFunction.Debt','LineWidth',lwidth)
axis tight;
set(gca,'FontSize',fsize)
title('Labor 2 Pol Fun')
xlabel('Debt','FontSize',fsize)

cd(path_Figure);
set(gcf,'PaperOrientation','landscape','PaperPosition',[0 0 16 12])
set(gcf, 'PaperSize', [16 12]);
print(figure5,'-dpdf',name5)
cd(path_Main);


figure6 = figure(6);
name6 = strcat('TangebleCapitalOnePolicyFunction_Chebychev.pdf');

subplot(2,2,1)
plot(vGridTangebleCapital',Q1.functions.vTangebleCapitalOnePolicyFunction.TangebleCapital','LineWidth',lwidth)
axis tight;
legend('both high','z1 high','z2 high','both low')
set(gca,'FontSize',fsize)
title('TangebleCapital 1 Pol Fun')
xlabel('Tangeble Capital','FontSize',fsize)
set(legend,'FontSize',fsize1)
set(legend,'Position',[0.22,0.85,0.01,0.01])
set(legend,'Box','off')


subplot(2,2,2)
plot(vGridIntangebleCapital',Q1.functions.vTangebleCapitalOnePolicyFunction.IntangebleCapital','LineWidth',lwidth)
axis tight;
set(gca,'FontSize',fsize)
title('Tangeble Capital 1 Pol Fun')
xlabel('Intangeble Capital','FontSize',fsize)

subplot(2,2,3)
plot(vGridDebt',Q1.functions.vTangebleCapitalOnePolicyFunction.Debt','LineWidth',lwidth)
axis tight;
set(gca,'FontSize',fsize)
title('Tangeble Capital 1 Pol Fun')
xlabel('Debt','FontSize',fsize)

cd(path_Figure);
set(gcf,'PaperOrientation','landscape','PaperPosition',[0 0 16 12])
set(gcf, 'PaperSize', [16 12]);
print(figure6,'-dpdf',name6)
cd(path_Main);


figure7 = figure(7);
name7 = strcat('IntangebleCapitalOnePolicyFunction_Chebychev.pdf');

subplot(2,2,1)
plot(vGridTangebleCapital',Q1.functions.vIntangebleCapitalOnePolicyFunction.TangebleCapital','LineWidth',lwidth)
axis tight;
legend('both high','z1 high','z2 high','both low')
set(gca,'FontSize',fsize)
title('IntangebleCapital 1 Pol Fun')
xlabel('Tangeble Capital','FontSize',fsize)
set(legend,'FontSize',fsize1)
set(legend,'Position',[0.22,0.85,0.01,0.01])
set(legend,'Box','off')


subplot(2,2,2)
plot(vGridIntangebleCapital',Q1.functions.vIntangebleCapitalOnePolicyFunction.IntangebleCapital','LineWidth',lwidth)
axis tight;
set(gca,'FontSize',fsize)
title('Intangeble Capital 1 Pol Fun')
xlabel('Intangeble Capital','FontSize',fsize)

subplot(2,2,3)
plot(vGridDebt',Q1.functions.vIntangebleCapitalOnePolicyFunction.Debt','LineWidth',lwidth)
axis tight;
set(gca,'FontSize',fsize)
title('Intangeble Capital 1 Pol Fun')
xlabel('Debt','FontSize',fsize)

cd(path_Figure);
set(gcf,'PaperOrientation','landscape','PaperPosition',[0 0 16 12])
set(gcf, 'PaperSize', [16 12]);
print(figure7,'-dpdf',name7)
cd(path_Main);

figure8 = figure(8);
name8 = strcat('DebtPolicyFunction_Chebychev.pdf');

subplot(2,2,1)
plot(vGridTangebleCapital',Q1.functions.vDebtPolicyFunction.TangebleCapital','LineWidth',lwidth)
axis tight;
legend('both high','z1 high','z2 high','both low')
set(gca,'FontSize',fsize)
title('Debt Pol Fun')
xlabel('Tangeble Capital','FontSize',fsize)
set(legend,'FontSize',fsize1)
set(legend,'Position',[0.22,0.85,0.01,0.01])
set(legend,'Box','off')


subplot(2,2,2)
plot(vGridIntangebleCapital',Q1.functions.vDebtPolicyFunction.IntangebleCapital','LineWidth',lwidth)
axis tight;
set(gca,'FontSize',fsize)
title('Debt Pol Fun')
xlabel('Intangeble Capital','FontSize',fsize)

subplot(2,2,3)
plot(vGridDebt',Q1.functions.vDebtPolicyFunction.Debt','LineWidth',lwidth)
axis tight;
set(gca,'FontSize',fsize)
title('Debt Pol Fun')
xlabel('Debt','FontSize',fsize)

cd(path_Figure);
set(gcf,'PaperOrientation','landscape','PaperPosition',[0 0 16 12])
set(gcf, 'PaperSize', [16 12]);
print(figure8,'-dpdf',name8)
cd(path_Main);

%% Simulation

% create markov path
% assume the economy starts from state 1

nGridTime = 1000;
vGridPeriod = [1:1:nGridTime+1]';
mStateRealization = zeros(nGridTime+1,2);
mSTateRealization(1,:) = [vProductivityT(1),vProductivityI(1)];

for nTime = 2:nGridTime+1
    vStateSim = rand([1,2]);
    if vStateSim(1) <= 0.95 % same state
        mStateRealization(nTime,1) = mStateRealization(nTime-1,1);
    else % diff state
        if mStateRealization(nTime-1,1) == vProductivityT(1)
            mStateRealization(nTime,1) = vProductivityT(2);
        else
            mStateRealization(nTime,1) = vProductivityT(1);
        end
    end
    
    if vStateSim(2) <= 0.90 % same state
        mStateRealization(nTime,2) = mStateRealization(nTime-1,2);
    else % diff state
        if mStateRealization(nTime-1,2) == vProductivityI(1)
            mStateRealization(nTime,2) = vProductivityI(2);
        else
            mStateRealization(nTime,2) = vProductivityI(1);
        end
    end
end

% save for replication
save('Result.mat','mStateRealization','-append');
load('Result.mat');


% Assume current state variables are at the deterministic steady state
% level

vSimulatedPath.TangebleCapital = zeros(nGridTime+1,1);
vSimulatedPath.TangebleCapital(1) = kt_ss;
vSimulatedPath.IntangebleCapital = zeros(nGridTime+1,1);
vSimulatedPath.IntangebleCapital(1) = ki_ss;
vSimulatedPath.TangebleCapitalOne = zeros(nGridTime+1,1);
vSimulatedPath.TangebleCapitalOne(1) = kt1_ss;
vSimulatedPath.IntangebleCapitalOne = zeros(nGridTime+1,1);
vSimulatedPath.IntangebleCapitalOne(1) = ki1_ss;
vSimulatedPath.LaborOne = zeros(nGridTime+1,1);
vSimulatedPath.LaborOne(1) = l1_ss;
vSimulatedPath.LaborTwo = zeros(nGridTime+1,1);
vSimulatedPath.LaborTwo(1) = l2_ss;
vSimulatedPath.Debt = zeros(nGridTime+1,1);
vSimulatedPath.Debt(1) = b_ss;
vSimulatedPath.Value = zeros(nGridTime+1,1);
vSimulatedPath.Value(1) = J_ss;

for nTime = 2:nGridTime+1
    % scale down state variables
    nTangebleCapitalTodayScaledown = 2*(vSimulatedPath.TangebleCapital(nTime-1)-nTangebleCapitalMin)/(nTangebleCapitalMax-nTangebleCapitalMin)-1;
    nIntangebleCapitalTodayScaledown = 2*(vSimulatedPath.IntangebleCapital(nTime-1)-nIntangebleCapitalMin)/(nIntangebleCapitalMax-nIntangebleCapitalMin)-1;
    nDebtTodayScaledown = 2*(vSimulatedPath.Debt(nTime-1)-nDebtMin)/(nDebtMax-nDebtMin)-1;
    
    % Chebychev
    vChebychev_kt_sim = ones(1,nMaxOrder);
    vChebychev_ki_sim = ones(1,nMaxOrder);
    vChebychev_b_sim  = ones(1,nMaxOrder);
    vChebychev_kt_sim(2) = nTangebleCapitalTodayScaledown;
    vChebychev_ki_sim(2) = nIntangebleCapitalTodayScaledown;
    vChebychev_b_sim(2)  = nDebtTodayScaledown;
    for i1=3:nMaxOrder
        vChebychev_kt_sim(i1) = 2*vChebychev_kt_sim(2)*vChebychev_kt_sim(i1-1)-vChebychev_kt_sim(i1-2);
        vChebychev_ki_sim(i1) = 2*vChebychev_ki_sim(2)*vChebychev_ki_sim(i1-1)-vChebychev_ki_sim(i1-2);
        vChebychev_b_sim(i1)  = 2*vChebychev_b_sim(2)*vChebychev_b_sim(i1-1)-vChebychev_b_sim(i1-2);
    end
    vChebychevToday = [1; vChebychev_kt_sim(2:nMaxOrder)'; vChebychev_ki_sim(2:nMaxOrder)'; vChebychev_b_sim(2:nMaxOrder)'; % const and no cross terms (13)
                       vChebychev_kt_sim(2)*vChebychev_ki_sim(2); vChebychev_kt_sim(2)*vChebychev_b_sim(2); vChebychev_ki_sim(2)*vChebychev_b_sim(2); % cross 1*1 (3)
                       vChebychev_kt_sim(2)*vChebychev_ki_sim(3); vChebychev_kt_sim(2)*vChebychev_b_sim(3); vChebychev_ki_sim(2)*vChebychev_kt_sim(3); 
                       vChebychev_ki_sim(2)*vChebychev_b_sim(3);  vChebychev_b_sim(2)*vChebychev_kt_sim(3); vChebychev_b_sim(2)*vChebychev_ki_sim(3); % cross 1*2 (6)
                       vChebychev_kt_sim(3)*vChebychev_ki_sim(3); vChebychev_kt_sim(3)*vChebychev_b_sim(3); vChebychev_ki_sim(3)*vChebychev_b_sim(3)];% cross 2*2 (3)
    
    if mStateRealization(nTime,1) == log(1.05)
        nProductivityT = 1;
    else
        nProductivityT = 2;
    end
    if mStateRealization(nTime,2) == log(1.10)
        nProductivityI = 1;
    else
        nProductivityI = 2;
    end
        
    
    nProductivity = (nProductivityT-1)*nGridProductivityT+nProductivityI;
    
    
    vSimulatedPath.TangebleCapital(nTime) = vec(mTheta(nIndexkt,:,nProductivity))'*vChebychevToday;
    vSimulatedPath.IntangebleCapital(nTime) = vec(mTheta(nIndexki,:,nProductivity))'*vChebychevToday;
    vSimulatedPath.LaborOne(nTime) = vec(mTheta(nIndexl1,:,nProductivity))'*vChebychevToday;
    vSimulatedPath.LaborTwo(nTime) = vec(mTheta(nIndexl2,:,nProductivity))'*vChebychevToday;
    nTangebleCapitalTwo = vec(mTheta(nIndexkt2,:,nProductivity))'*vChebychevToday;
    nIntangebleCapitalTwo = vec(mTheta(nIndexki2,:,nProductivity))'*vChebychevToday;
    vSimulatedPath.TangebleCapitalOne(nTime) = vSimulatedPath.TangebleCapital(nTime-1)-nTangebleCapitalTwo;
    vSimulatedPath.IntangebleCapitalOne(nTime) = vSimulatedPath.IntangebleCapital(nTime-1)-nIntangebleCapitalTwo;
    vSimulatedPath.Debt(nTime) = vec(mTheta(nIndexb,:,nProductivity))'*vChebychevToday;
    vSimulatedPath.Value(nTime) = vec(mTheta(nIndexJ,:,nProductivity))'*vChebychevToday;
end
    
    
name99 = strcat('Dynamics_state_Chebychev.pdf');
figure99 = figure(99);
subplot(3,2,1)
plot(vGridPeriod,vSimulatedPath.TangebleCapital)
title('Tangeble Capital')
axis tight;
set(gca,'FontSize',fsize)

subplot(3,2,2)
plot(vGridPeriod,vSimulatedPath.IntangebleCapital)
title('Intangeble Capital')
axis tight;
set(gca,'FontSize',fsize)

subplot(3,2,3)
plot(vGridPeriod,vSimulatedPath.Debt)
title('Debt')
axis tight;
set(gca,'FontSize',fsize)

subplot(3,2,4)
plot(vGridPeriod,mStateRealization(:,1))
title('z1')
axis tight;
set(gca,'FontSize',fsize)

subplot(3,2,5)
plot(vGridPeriod,mStateRealization(:,2))
title('z2')
axis tight;
set(gca,'FontSize',fsize)

cd(path_Figure);
set(gcf,'PaperOrientation','landscape','PaperPosition',[0 0 15 10])
set(gcf, 'PaperSize', [15 10]);
print(figure99,'-dpdf',name99)
cd(path_Main);


name100 = strcat('Dynamics_endo_Chebychev.pdf');
figure100 = figure(100);
subplot(3,2,1)
plot(vGridPeriod,vSimulatedPath.LaborOne)
title('Labor 1')
axis tight;
set(gca,'FontSize',fsize)

subplot(3,2,2)
plot(vGridPeriod,vSimulatedPath.LaborTwo)
title('Labor 2')
axis tight;
set(gca,'FontSize',fsize)

subplot(3,2,3)
plot(vGridPeriod,vSimulatedPath.TangebleCapitalOne)
title('Tangeble Capital 1')
axis tight;
set(gca,'FontSize',fsize)

subplot(3,2,4)
plot(vGridPeriod,vSimulatedPath.IntangebleCapital)
title('Intangeble Capital 1')
axis tight;
set(gca,'FontSize',fsize)

subplot(3,2,5)
plot(vGridPeriod,vSimulatedPath.IntangebleCapital)
title('Value Function')
axis tight;
set(gca,'FontSize',fsize)





cd(path_Figure);
set(gcf,'PaperOrientation','landscape','PaperPosition',[0 0 15 10])
set(gcf, 'PaperSize', [15 10]);
print(figure100,'-dpdf',name100)
cd(path_Main);





