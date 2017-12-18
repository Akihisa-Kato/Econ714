% Code for Econ714 HW2 Q2
% Firm's decision


%% Preliminary
%clear all;
close all; clc;

fsize = 7;
fsize1 = 5;
lwidth=1.5;

path_Main=pwd;
cd ../.
path_Figure=[pwd '/Main/Figure'];

addpath(path_Main);
addpath(genpath(path_Figure));

cd(path_Main);


resultmat = strcat('Result.mat');

%% Run Dynare
dynare Q2.mod;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  1. Value Function 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Tangeble Capital

% find index that value function is placed in y vector
nIndexValueFunction = find(strcmp('J', outDynare.label_y));
nIndexTangebleCapitalOne = find(strcmp('kt1', outDynare.label_y));
nIndexIntangebleCapitalOne = find(strcmp('ki1', outDynare.label_y));
nIndexLaborOne = find(strcmp('l1', outDynare.label_y));
nIndexLaborTwo = find(strcmp('l2', outDynare.label_y));


nIndexTangebleCapital = find(strcmp('kt', outDynare.label_v));
nIndexIntangebleCapital = find(strcmp('ki', outDynare.label_v));
nIndexZ1 = find(strcmp('z1', outDynare.label_v));
nIndexZ2 = find(strcmp('z2', outDynare.label_v));
nIndexez1 = find(strcmp('ez1', outDynare.label_v));
nIndexez2 = find(strcmp('ez2', outDynare.label_v));
nIndexDebt = find(strcmp('b', outDynare.label_v));



% define coefficients for Control
vSteadyStateControl = outDynare.gSteadyState;
nControlVariable = size(vSteadyStateControl);
vFirstCoeffStateOnControl = outDynare.gv;
mSecondCoeffStateOnControl = outDynare.gvv/2;
mThridCoeffStateOnControl = outDynare.gvvv/6;
mSecondCoeffPerturbationOnControl = outDynare.gss/2;
mThirdCoeffCrossTermOnControl = outDynare.gssv/6;
mThirdCoeffPerturbationOnControl = outDynare.gsss/6; % Theoretially 0's

% define coefficients for State
vFirstCoeffStateOnState = outDynare.hv;
mSecondCoeffStateOnState = outDynare.hvv/2;
mThridCoeffStateOnState = outDynare.hvvv/6;
mSecondCoeffPerturbationOnState = outDynare.hss/2;
mThirdCoeffCrossTermOnState = outDynare.hssv/6;
mThirdCoeffPerturbationOnState = outDynare.hsss/6; % Theoretially 0's



% Define Deviations from Steady State
vSteadyStateState = outDynare.hSteadyState;
nStateVariable = size(vSteadyStateState);


nMaxTangebleCapital = 1.3*vSteadyStateState(nIndexTangebleCapital);
nMinTangebleCapital = 0.7*vSteadyStateState(nIndexTangebleCapital);
nGridTangebleCapital = 150;
vGridTangebleCapital = linspace(nMinTangebleCapital,nMaxTangebleCapital,nGridTangebleCapital);


% No shock
for nTangebleCapital = 1:nGridTangebleCapital
    vStateVariable = zeros(nStateVariable);
    vStateVariable(nIndexTangebleCapital) = vGridTangebleCapital(nTangebleCapital)-vSteadyStateState(nIndexTangebleCapital);
    vValueFunction.TangebleCapital(5,nTangebleCapital) = vSteadyStateControl(nIndexValueFunction)+vFirstCoeffStateOnControl(nIndexValueFunction,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexValueFunction)...
        +mSecondCoeffStateOnControl(nIndexValueFunction,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexValueFunction,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexValueFunction,:)*vStateVariable;
end



nGridShockZ1 = 2;
nGridShockZ2 = 2;

vGridShockZ1 = [1.05 0.95];
vGridShockZ2 = [1.10 0.90];

for nShockZ2 = 1:nGridShockZ1
    for nShockZ1 = 1:nGridShockZ2
        nShock = (nShockZ1-1)*nGridShockZ2+nShockZ2;
        vStateVariable = zeros(nStateVariable);
        vStateVariable(nIndexZ1) = log(vGridShockZ1(nShockZ1));
        vStateVariable(nIndexZ2) = log(vGridShockZ2(nShockZ2));
                
        % both high
        for nTangebleCapital = 1:nGridTangebleCapital
            vStateVariable(nIndexTangebleCapital) = vGridTangebleCapital(nTangebleCapital)-vSteadyStateState(nIndexTangebleCapital);
            vValueFunction.TangebleCapital(nShock,nTangebleCapital) = vSteadyStateControl(nIndexValueFunction)+vFirstCoeffStateOnControl(nIndexValueFunction,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexValueFunction)...
                +mSecondCoeffStateOnControl(nIndexValueFunction,:)*kron(vStateVariable,vStateVariable)+vec(mThridCoeffStateOnControl(nIndexValueFunction,:,:,:))'*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexValueFunction,:)*vStateVariable;
        end
    end
end
%% Intangeble Capital


nMaxIntangebleCapital = 1.3*vSteadyStateState(nIndexIntangebleCapital);
nMinIntangebleCapital = 0.7*vSteadyStateState(nIndexIntangebleCapital);
nGridIntangebleCapital = 150;
vGridIntangebleCapital = linspace(nMinIntangebleCapital,nMaxIntangebleCapital,nGridIntangebleCapital);


% No shock
for nIntangebleCapital = 1:nGridIntangebleCapital
    vStateVariable = zeros(nStateVariable);
    vStateVariable(nIndexIntangebleCapital) = vGridIntangebleCapital(nIntangebleCapital)-vSteadyStateState(nIndexIntangebleCapital);
    vValueFunction.IntangebleCapital(5,nIntangebleCapital) = vSteadyStateControl(nIndexValueFunction)+vFirstCoeffStateOnControl(nIndexValueFunction,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexValueFunction)...
        +mSecondCoeffStateOnControl(nIndexValueFunction,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexValueFunction,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexValueFunction,:)*vStateVariable;
end


for nShockZ2 = 1:nGridShockZ1
    for nShockZ1 = 1:nGridShockZ2
        nShock = (nShockZ1-1)*nGridShockZ2+nShockZ2;
        vStateVariable = zeros(nStateVariable);
        vStateVariable(nIndexZ1) = log(vGridShockZ1(nShockZ1));
        vStateVariable(nIndexZ2) = log(vGridShockZ2(nShockZ2));
        
        for nIntangebleCapital = 1:nGridIntangebleCapital
            vStateVariable(nIndexIntangebleCapital) = vGridIntangebleCapital(nIntangebleCapital)-vSteadyStateState(nIndexIntangebleCapital);
            vValueFunction.IntangebleCapital(nShock,nIntangebleCapital) = vSteadyStateControl(nIndexValueFunction)+vFirstCoeffStateOnControl(nIndexValueFunction,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexValueFunction)...
                +mSecondCoeffStateOnControl(nIndexValueFunction,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexValueFunction,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexValueFunction,:)*vStateVariable;
        end
    end
end

%% Debt


nMaxDebt = 0.7*vSteadyStateState(nIndexDebt);
nMinDebt = 1.3*vSteadyStateState(nIndexDebt);
nGridDebt = 150;
vGridDebt = linspace(nMinDebt,nMaxDebt,nGridDebt);


% No shock
for nDebt = 1:nGridDebt
    vStateVariable = zeros(nStateVariable);
    vStateVariable(nIndexDebt) = vGridDebt(nDebt)-vSteadyStateState(nIndexDebt);
    vValueFunction.Debt(5,nDebt) = vSteadyStateControl(nIndexValueFunction)+vFirstCoeffStateOnControl(nIndexValueFunction,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexValueFunction)...
        +mSecondCoeffStateOnControl(nIndexValueFunction,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexValueFunction,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexValueFunction,:)*vStateVariable;
end


for nShockZ2 = 1:nGridShockZ1
    for nShockZ1 = 1:nGridShockZ2
        nShock = (nShockZ1-1)*nGridShockZ2+nShockZ2;
        vStateVariable = zeros(nStateVariable);
        vStateVariable(nIndexZ1) = log(vGridShockZ1(nShockZ1));
        vStateVariable(nIndexZ2) = log(vGridShockZ2(nShockZ2));
        
        for nDebt = 1:nGridDebt
            vStateVariable(nIndexDebt) = vGridDebt(nDebt)-vSteadyStateState(nIndexDebt);
            vValueFunction.Debt(nShock,nDebt) = vSteadyStateControl(nIndexValueFunction)+vFirstCoeffStateOnControl(nIndexValueFunction,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexValueFunction)...
                +mSecondCoeffStateOnControl(nIndexValueFunction,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexValueFunction,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexValueFunction,:)*vStateVariable;
        end
    end
end

save(resultmat,'vValueFunction','-append')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  2. Policy Functoins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%-----------------------------------
% Tangeble Capital Policy Function 
%-----------------------------------
% Changing Tangeble capital

% No shock
for nTangebleCapital = 1:nGridTangebleCapital
    vStateVariable = zeros(nStateVariable);
    vStateVariable(nIndexTangebleCapital) = vGridTangebleCapital(nTangebleCapital)-vSteadyStateState(nIndexTangebleCapital);
    vPolicyTangebleCapital.TangebleCapital(5,nTangebleCapital) = vSteadyStateState(nIndexTangebleCapital)+vFirstCoeffStateOnState(nIndexTangebleCapital,:)*vStateVariable+mSecondCoeffPerturbationOnState(nIndexTangebleCapital)...
        +mSecondCoeffStateOnState(nIndexTangebleCapital,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnState(nIndexTangebleCapital,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnState(nIndexTangebleCapital,:)*vStateVariable;
end

for nShockZ2 = 1:nGridShockZ1
    for nShockZ1 = 1:nGridShockZ2
        nShock = (nShockZ1-1)*nGridShockZ2+nShockZ2;
        vStateVariable = zeros(nStateVariable);
        vStateVariable(nIndexZ1) = log(vGridShockZ1(nShockZ1));
        vStateVariable(nIndexZ2) = log(vGridShockZ2(nShockZ2));
        
        for nTangebleCapital = 1:nGridTangebleCapital
            vStateVariable(nIndexTangebleCapital) = vGridTangebleCapital(nTangebleCapital)-vSteadyStateState(nIndexTangebleCapital);
            vPolicyTangebleCapital.TangebleCapital(nShock,nTangebleCapital) = vSteadyStateState(nIndexTangebleCapital)+vFirstCoeffStateOnState(nIndexTangebleCapital,:)*vStateVariable+mSecondCoeffPerturbationOnState(nIndexTangebleCapital)...
                +mSecondCoeffStateOnState(nIndexTangebleCapital,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnState(nIndexTangebleCapital,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnState(nIndexTangebleCapital,:)*vStateVariable;
        end
    end
end


% Changing Intangeble capital

% No shock
for nIntangebleCapital = 1:nGridIntangebleCapital
    vStateVariable = zeros(nStateVariable);
    vStateVariable(nIndexIntangebleCapital) = vGridIntangebleCapital(nIntangebleCapital)-vSteadyStateState(nIndexIntangebleCapital);
    vPolicyTangebleCapital.IntangebleCapital(5,nIntangebleCapital) = vSteadyStateState(nIndexTangebleCapital)+vFirstCoeffStateOnState(nIndexTangebleCapital,:)*vStateVariable+mSecondCoeffPerturbationOnState(nIndexTangebleCapital)...
        +mSecondCoeffStateOnState(nIndexTangebleCapital,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnState(nIndexTangebleCapital,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnState(nIndexTangebleCapital,:)*vStateVariable;
end

for nShockZ2 = 1:nGridShockZ1
    for nShockZ1 = 1:nGridShockZ2
        nShock = (nShockZ1-1)*nGridShockZ2+nShockZ2;
        vStateVariable = zeros(nStateVariable);
        vStateVariable(nIndexZ1) = log(vGridShockZ1(nShockZ1));
        vStateVariable(nIndexZ2) = log(vGridShockZ2(nShockZ2));
        for nIntangebleCapital = 1:nGridIntangebleCapital
            vStateVariable(nIndexIntangebleCapital) = vGridIntangebleCapital(nIntangebleCapital)-vSteadyStateState(nIndexIntangebleCapital);
            vPolicyTangebleCapital.IntangebleCapital(nShock,nIntangebleCapital) = vSteadyStateState(nIndexTangebleCapital)+vFirstCoeffStateOnState(nIndexTangebleCapital,:)*vStateVariable+mSecondCoeffPerturbationOnState(nIndexTangebleCapital)...
                +mSecondCoeffStateOnState(nIndexTangebleCapital,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnState(nIndexTangebleCapital,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnState(nIndexTangebleCapital,:)*vStateVariable;
        end
    end
end


% Changing Debt

% No shock
for nDebt = 1:nGridDebt
    vStateVariable = zeros(nStateVariable);
    vStateVariable(nIndexDebt) = vGridDebt(nDebt)-vSteadyStateState(nIndexDebt);
    vPolicyTangebleCapital.Debt(5,nDebt) = vSteadyStateState(nIndexTangebleCapital)+vFirstCoeffStateOnState(nIndexTangebleCapital,:)*vStateVariable+mSecondCoeffPerturbationOnState(nIndexTangebleCapital)...
        +mSecondCoeffStateOnState(nIndexTangebleCapital,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnState(nIndexTangebleCapital,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnState(nIndexTangebleCapital,:)*vStateVariable;
end

for nShockZ2 = 1:nGridShockZ1
    for nShockZ1 = 1:nGridShockZ2
        nShock = (nShockZ1-1)*nGridShockZ2+nShockZ2;
        vStateVariable = zeros(nStateVariable);
        vStateVariable(nIndexZ1) = log(vGridShockZ1(nShockZ1));
        vStateVariable(nIndexZ2) = log(vGridShockZ2(nShockZ2));
        for nDebt = 1:nGridDebt
            vStateVariable(nIndexDebt) = vGridDebt(nDebt)-vSteadyStateState(nIndexDebt);
            vPolicyTangebleCapital.Debt(nShock,nDebt) = vSteadyStateState(nIndexTangebleCapital)+vFirstCoeffStateOnState(nIndexTangebleCapital,:)*vStateVariable+mSecondCoeffPerturbationOnState(nIndexTangebleCapital)...
                +mSecondCoeffStateOnState(nIndexTangebleCapital,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnState(nIndexTangebleCapital,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnState(nIndexTangebleCapital,:)*vStateVariable;
        end
    end
end

save(resultmat,'vPolicyTangebleCapital','-append')

%%
%-------------------------------------- 
% Inangeble Capital Policy Function 
%---------------------------------------

% Changing Tangeble capital

% No shock
for nTangebleCapital = 1:nGridTangebleCapital
    vStateVariable = zeros(nStateVariable);
    vStateVariable(nIndexTangebleCapital) = vGridTangebleCapital(nTangebleCapital)-vSteadyStateState(nIndexTangebleCapital);
    vPolicyIntangebleCapital.TangebleCapital(5,nTangebleCapital) = vSteadyStateState(nIndexIntangebleCapital)+vFirstCoeffStateOnState(nIndexIntangebleCapital,:)*vStateVariable+mSecondCoeffPerturbationOnState(nIndexIntangebleCapital)...
        +mSecondCoeffStateOnState(nIndexIntangebleCapital,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnState(nIndexIntangebleCapital,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnState(nIndexIntangebleCapital,:)*vStateVariable;
end

for nShockZ2 = 1:nGridShockZ1
    for nShockZ1 = 1:nGridShockZ2
        nShock = (nShockZ1-1)*nGridShockZ2+nShockZ2;
        vStateVariable = zeros(nStateVariable);
        vStateVariable(nIndexZ1) = log(vGridShockZ1(nShockZ1));
        vStateVariable(nIndexZ2) = log(vGridShockZ2(nShockZ2));
        
        for nTangebleCapital = 1:nGridTangebleCapital
            vStateVariable(nIndexTangebleCapital) = vGridTangebleCapital(nTangebleCapital)-vSteadyStateState(nIndexTangebleCapital);
            vPolicyIntangebleCapital.TangebleCapital(nShock,nTangebleCapital) = vSteadyStateState(nIndexIntangebleCapital)+vFirstCoeffStateOnState(nIndexIntangebleCapital,:)*vStateVariable+mSecondCoeffPerturbationOnState(nIndexIntangebleCapital)...
                +mSecondCoeffStateOnState(nIndexIntangebleCapital,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnState(nIndexIntangebleCapital,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnState(nIndexIntangebleCapital,:)*vStateVariable;
        end
    end
end


% Changing Intangeble capital

% No shock
for nIntangebleCapital = 1:nGridIntangebleCapital
    vStateVariable = zeros(nStateVariable);
    vStateVariable(nIndexIntangebleCapital) = vGridIntangebleCapital(nIntangebleCapital)-vSteadyStateState(nIndexIntangebleCapital);
    vPolicyIntangebleCapital.IntangebleCapital(5,nIntangebleCapital) = vSteadyStateState(nIndexIntangebleCapital)+vFirstCoeffStateOnState(nIndexIntangebleCapital,:)*vStateVariable+mSecondCoeffPerturbationOnState(nIndexIntangebleCapital)...
        +mSecondCoeffStateOnState(nIndexIntangebleCapital,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnState(nIndexIntangebleCapital,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnState(nIndexIntangebleCapital,:)*vStateVariable;
end

for nShockZ2 = 1:nGridShockZ1
    for nShockZ1 = 1:nGridShockZ2
        nShock = (nShockZ1-1)*nGridShockZ2+nShockZ2;
        vStateVariable = zeros(nStateVariable);
        vStateVariable(nIndexZ1) = log(vGridShockZ1(nShockZ1));
        vStateVariable(nIndexZ2) = log(vGridShockZ2(nShockZ2));
        for nIntangebleCapital = 1:nGridIntangebleCapital
            vStateVariable(nIndexIntangebleCapital) = vGridIntangebleCapital(nIntangebleCapital)-vSteadyStateState(nIndexIntangebleCapital);
            vPolicyIntangebleCapital.IntangebleCapital(nShock,nIntangebleCapital) = vSteadyStateState(nIndexIntangebleCapital)+vFirstCoeffStateOnState(nIndexIntangebleCapital,:)*vStateVariable+mSecondCoeffPerturbationOnState(nIndexIntangebleCapital)...
                +mSecondCoeffStateOnState(nIndexIntangebleCapital,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnState(nIndexIntangebleCapital,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnState(nIndexIntangebleCapital,:)*vStateVariable;
        end
    end
end

% Changing Debt
% No shock
for nDebt = 1:nGridDebt
    vStateVariable = zeros(nStateVariable);
    vStateVariable(nIndexDebt) = vGridDebt(nDebt)-vSteadyStateState(nIndexDebt);
    vPolicyIntangebleCapital.Debt(5,nDebt) = vSteadyStateState(nIndexIntangebleCapital)+vFirstCoeffStateOnState(nIndexIntangebleCapital,:)*vStateVariable+mSecondCoeffPerturbationOnState(nIndexIntangebleCapital)...
        +mSecondCoeffStateOnState(nIndexIntangebleCapital,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnState(nIndexIntangebleCapital,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnState(nIndexIntangebleCapital,:)*vStateVariable;
end

for nShockZ2 = 1:nGridShockZ1
    for nShockZ1 = 1:nGridShockZ2
        nShock = (nShockZ1-1)*nGridShockZ2+nShockZ2;
        vStateVariable = zeros(nStateVariable);
        vStateVariable(nIndexZ1) = log(vGridShockZ1(nShockZ1));
        vStateVariable(nIndexZ2) = log(vGridShockZ2(nShockZ2));
        
        for nDebt = 1:nGridDebt
            vStateVariable(nIndexDebt) = vGridDebt(nDebt)-vSteadyStateState(nIndexDebt);
            vPolicyIntangebleCapital.Debt(nShock,nDebt) = vSteadyStateState(nIndexIntangebleCapital)+vFirstCoeffStateOnState(nIndexIntangebleCapital,:)*vStateVariable+mSecondCoeffPerturbationOnState(nIndexIntangebleCapital)...
                +mSecondCoeffStateOnState(nIndexIntangebleCapital,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnState(nIndexIntangebleCapital,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnState(nIndexIntangebleCapital,:)*vStateVariable;     

        end
    end
end

save(resultmat,'vPolicyIntangebleCapital','-append')

%%
%-------------------------------------- 
% Debt Policy Function 
%---------------------------------------

% Changing Tangeble capital

% No shock
for nTangebleCapital = 1:nGridTangebleCapital
    vStateVariable = zeros(nStateVariable);
    vStateVariable(nIndexTangebleCapital) = vGridTangebleCapital(nTangebleCapital)-vSteadyStateState(nIndexTangebleCapital);
    vPolicyDebt.TangebleCapital(5,nTangebleCapital) = vSteadyStateState(nIndexDebt)+vFirstCoeffStateOnState(nIndexDebt,:)*vStateVariable+mSecondCoeffPerturbationOnState(nIndexDebt)...
        +mSecondCoeffStateOnState(nIndexDebt,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnState(nIndexDebt,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnState(nIndexDebt,:)*vStateVariable;
end

for nShockZ2 = 1:nGridShockZ1
    for nShockZ1 = 1:nGridShockZ2
        nShock = (nShockZ1-1)*nGridShockZ2+nShockZ2;
        vStateVariable = zeros(nStateVariable);
        vStateVariable(nIndexZ1) = log(vGridShockZ1(nShockZ1));
        vStateVariable(nIndexZ2) = log(vGridShockZ2(nShockZ2));
        
        for nTangebleCapital = 1:nGridTangebleCapital
            vStateVariable(nIndexTangebleCapital) = vGridTangebleCapital(nTangebleCapital)-vSteadyStateState(nIndexTangebleCapital);
            vPolicyDebt.TangebleCapital(nShock,nTangebleCapital) = vSteadyStateState(nIndexDebt)+vFirstCoeffStateOnState(nIndexDebt,:)*vStateVariable+mSecondCoeffPerturbationOnState(nIndexDebt)...
                +mSecondCoeffStateOnState(nIndexDebt,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnState(nIndexDebt,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnState(nIndexDebt,:)*vStateVariable;
        end
    end
end


% Changing Intangeble capital

% No shock
for nIntangebleCapital = 1:nGridIntangebleCapital
    vStateVariable = zeros(nStateVariable);
    vStateVariable(nIndexIntangebleCapital) = vGridIntangebleCapital(nIntangebleCapital)-vSteadyStateState(nIndexIntangebleCapital);
    vPolicyDebt.IntangebleCapital(5,nIntangebleCapital) = vSteadyStateState(nIndexDebt)+vFirstCoeffStateOnState(nIndexDebt,:)*vStateVariable+mSecondCoeffPerturbationOnState(nIndexDebt)...
        +mSecondCoeffStateOnState(nIndexDebt,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnState(nIndexDebt,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnState(nIndexDebt,:)*vStateVariable;
end

for nShockZ2 = 1:nGridShockZ1
    for nShockZ1 = 1:nGridShockZ2
        nShock = (nShockZ1-1)*nGridShockZ2+nShockZ2;
        vStateVariable = zeros(nStateVariable);
        vStateVariable(nIndexZ1) = log(vGridShockZ1(nShockZ1));
        vStateVariable(nIndexZ2) = log(vGridShockZ2(nShockZ2));
        for nIntangebleCapital = 1:nGridIntangebleCapital
            vStateVariable(nIndexIntangebleCapital) = vGridIntangebleCapital(nIntangebleCapital)-vSteadyStateState(nIndexIntangebleCapital);
            vPolicyDebt.IntangebleCapital(nShock,nIntangebleCapital) = vSteadyStateState(nIndexDebt)+vFirstCoeffStateOnState(nIndexDebt,:)*vStateVariable+mSecondCoeffPerturbationOnState(nIndexDebt)...
                +mSecondCoeffStateOnState(nIndexDebt,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnState(nIndexDebt,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnState(nIndexDebt,:)*vStateVariable;
        end
    end
end

% Changing Debt
% No shock
for nDebt = 1:nGridDebt
    vStateVariable = zeros(nStateVariable);
    vStateVariable(nIndexDebt) = vGridDebt(nDebt)-vSteadyStateState(nIndexDebt);
    vPolicyDebt.Debt(5,nDebt) = vSteadyStateState(nIndexDebt)+vFirstCoeffStateOnState(nIndexDebt,:)*vStateVariable+mSecondCoeffPerturbationOnState(nIndexDebt)...
        +mSecondCoeffStateOnState(nIndexDebt,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnState(nIndexDebt,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnState(nIndexDebt,:)*vStateVariable;
end

for nShockZ2 = 1:nGridShockZ1
    for nShockZ1 = 1:nGridShockZ2
        nShock = (nShockZ1-1)*nGridShockZ2+nShockZ2;
        vStateVariable = zeros(nStateVariable);
        vStateVariable(nIndexZ1) = log(vGridShockZ1(nShockZ1));
        vStateVariable(nIndexZ2) = log(vGridShockZ2(nShockZ2));
        
        for nDebt = 1:nGridDebt
            vStateVariable(nIndexDebt) = vGridDebt(nDebt)-vSteadyStateState(nIndexDebt);
            vPolicyDebt.Debt(nShock,nDebt) = vSteadyStateState(nIndexDebt)+vFirstCoeffStateOnState(nIndexDebt,:)*vStateVariable+mSecondCoeffPerturbationOnState(nIndexDebt)...
                +mSecondCoeffStateOnState(nIndexDebt,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnState(nIndexDebt,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnState(nIndexDebt,:)*vStateVariable;

        end
    end
end

save(resultmat,'vPolicyIntangebleCapital','-append')



%%
%-------------------------------------- 
% Tangeble Capital 1 Policy Function 
%---------------------------------------

% Changing Tangeble capital

% No shock
for nTangebleCapital = 1:nGridTangebleCapital
    vStateVariable = zeros(nStateVariable);
    vStateVariable(nIndexTangebleCapital) = vGridTangebleCapital(nTangebleCapital)-vSteadyStateState(nIndexTangebleCapital);
    vPolicyTangebleCapitalOne.TangebleCapital(5,nTangebleCapital) = vSteadyStateControl(nIndexTangebleCapitalOne)+vFirstCoeffStateOnControl(nIndexTangebleCapitalOne,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexTangebleCapitalOne)...
        +mSecondCoeffStateOnControl(nIndexTangebleCapitalOne,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexTangebleCapitalOne,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexTangebleCapitalOne,:)*vStateVariable;
end

for nShockZ2 = 1:nGridShockZ1
    for nShockZ1 = 1:nGridShockZ2
        nShock = (nShockZ1-1)*nGridShockZ2+nShockZ2;
        vStateVariable = zeros(nStateVariable);
        vStateVariable(nIndexZ1) = log(vGridShockZ1(nShockZ1));
        vStateVariable(nIndexZ2) = log(vGridShockZ2(nShockZ2));
        
        for nTangebleCapital = 1:nGridTangebleCapital
            vStateVariable(nIndexTangebleCapital) = vGridTangebleCapital(nTangebleCapital)-vSteadyStateState(nIndexTangebleCapital);
            vPolicyTangebleCapitalOne.TangebleCapital(nShock,nTangebleCapital) = vSteadyStateControl(nIndexTangebleCapitalOne)+vFirstCoeffStateOnControl(nIndexTangebleCapitalOne,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexTangebleCapitalOne)...
                +mSecondCoeffStateOnControl(nIndexTangebleCapitalOne,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexTangebleCapitalOne,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexTangebleCapitalOne,:)*vStateVariable;
        end
    end
end


% Changing Intangeble capital

% No shock
for nIntangebleCapital = 1:nGridIntangebleCapital
    vStateVariable = zeros(nStateVariable);
    vStateVariable(nIndexIntangebleCapital) = vGridIntangebleCapital(nIntangebleCapital)-vSteadyStateState(nIndexIntangebleCapital);
    vPolicyTangebleCapitalOne.IntangebleCapital(5,nIntangebleCapital) = vSteadyStateControl(nIndexTangebleCapitalOne)+vFirstCoeffStateOnControl(nIndexTangebleCapitalOne,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexTangebleCapitalOne)...
        +mSecondCoeffStateOnControl(nIndexTangebleCapitalOne,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexTangebleCapitalOne,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexTangebleCapitalOne,:)*vStateVariable;
end

for nShockZ2 = 1:nGridShockZ1
    for nShockZ1 = 1:nGridShockZ2
        nShock = (nShockZ1-1)*nGridShockZ2+nShockZ2;
        vStateVariable = zeros(nStateVariable);
        vStateVariable(nIndexZ1) = log(vGridShockZ1(nShockZ1));
        vStateVariable(nIndexZ2) = log(vGridShockZ2(nShockZ2));
        for nIntangebleCapital = 1:nGridIntangebleCapital
            vStateVariable(nIndexIntangebleCapital) = vGridIntangebleCapital(nIntangebleCapital)-vSteadyStateState(nIndexIntangebleCapital);
            vPolicyTangebleCapitalOne.IntangebleCapital(nShock,nIntangebleCapital) = vSteadyStateControl(nIndexTangebleCapitalOne)+vFirstCoeffStateOnControl(nIndexTangebleCapitalOne,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexTangebleCapitalOne)...
                +mSecondCoeffStateOnControl(nIndexTangebleCapitalOne,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexTangebleCapitalOne,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexTangebleCapitalOne,:)*vStateVariable;
        end
    end
end


% Changing Debt

% No shock
for nDebt = 1:nGridDebt
    vStateVariable = zeros(nStateVariable);
    vStateVariable(nIndexDebt) = vGridDebt(nDebt)-vSteadyStateState(nIndexDebt);
    vPolicyTangebleCapitalOne.Debt(5,nDebt) = vSteadyStateControl(nIndexTangebleCapitalOne)+vFirstCoeffStateOnControl(nIndexTangebleCapitalOne,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexTangebleCapitalOne)...
        +mSecondCoeffStateOnControl(nIndexTangebleCapitalOne,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexTangebleCapitalOne,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexTangebleCapitalOne,:)*vStateVariable;
end

for nShockZ2 = 1:nGridShockZ1
    for nShockZ1 = 1:nGridShockZ2
        nShock = (nShockZ1-1)*nGridShockZ2+nShockZ2;
        vStateVariable = zeros(nStateVariable);
        vStateVariable(nIndexZ1) = log(vGridShockZ1(nShockZ1));
        vStateVariable(nIndexZ2) = log(vGridShockZ2(nShockZ2));
        
        for nDebt = 1:nGridDebt
            vStateVariable(nIndexDebt) = vGridDebt(nDebt)-vSteadyStateState(nIndexDebt);
            vPolicyTangebleCapitalOne.Debt(nShock,nDebt) = vSteadyStateControl(nIndexTangebleCapitalOne)+vFirstCoeffStateOnControl(nIndexTangebleCapitalOne,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexTangebleCapitalOne)...
                +mSecondCoeffStateOnControl(nIndexTangebleCapitalOne,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexTangebleCapitalOne,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexTangebleCapitalOne,:)*vStateVariable;
        end
    end
end


save(resultmat,'vPolicyTangebleCapitalOne','-append')


%%
%-------------------------------------- 
% Intangeble Capital 1 Policy Function 
%---------------------------------------

% Changing Tangeble capital

% No shock
for nTangebleCapital = 1:nGridTangebleCapital
    vStateVariable = zeros(nStateVariable);
    vStateVariable(nIndexTangebleCapital) = vGridTangebleCapital(nTangebleCapital)-vSteadyStateState(nIndexTangebleCapital);
    vPolicyIntangebleCapitalOne.TangebleCapital(5,nTangebleCapital) = vSteadyStateControl(nIndexIntangebleCapitalOne)+vFirstCoeffStateOnControl(nIndexIntangebleCapitalOne,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexIntangebleCapitalOne)...
        +mSecondCoeffStateOnControl(nIndexIntangebleCapitalOne,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexIntangebleCapitalOne,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexIntangebleCapitalOne,:)*vStateVariable;
end

for nShockZ2 = 1:nGridShockZ1
    for nShockZ1 = 1:nGridShockZ2
        nShock = (nShockZ1-1)*nGridShockZ2+nShockZ2;
        vStateVariable = zeros(nStateVariable);
        vStateVariable(nIndexZ1) = log(vGridShockZ1(nShockZ1));
        vStateVariable(nIndexZ2) = log(vGridShockZ2(nShockZ2));
        
        for nTangebleCapital = 1:nGridTangebleCapital
            vStateVariable(nIndexTangebleCapital) = vGridTangebleCapital(nTangebleCapital)-vSteadyStateState(nIndexTangebleCapital);
            vPolicyIntangebleCapitalOne.TangebleCapital(nShock,nTangebleCapital) = vSteadyStateControl(nIndexIntangebleCapitalOne)+vFirstCoeffStateOnControl(nIndexIntangebleCapitalOne,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexIntangebleCapitalOne)...
                +mSecondCoeffStateOnControl(nIndexIntangebleCapitalOne,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexIntangebleCapitalOne,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexIntangebleCapitalOne,:)*vStateVariable;
        end
    end
end


% Changing Intangeble capital

% No shock
for nIntangebleCapital = 1:nGridIntangebleCapital
    vStateVariable = zeros(nStateVariable);
    vStateVariable(nIndexIntangebleCapital) = vGridIntangebleCapital(nIntangebleCapital)-vSteadyStateState(nIndexIntangebleCapital);
    vPolicyIntangebleCapitalOne.IntangebleCapital(5,nIntangebleCapital) = vSteadyStateControl(nIndexIntangebleCapitalOne)+vFirstCoeffStateOnControl(nIndexIntangebleCapitalOne,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexIntangebleCapitalOne)...
        +mSecondCoeffStateOnControl(nIndexIntangebleCapitalOne,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexIntangebleCapitalOne,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexIntangebleCapitalOne,:)*vStateVariable;
end

for nShockZ2 = 1:nGridShockZ1
    for nShockZ1 = 1:nGridShockZ2
        nShock = (nShockZ1-1)*nGridShockZ2+nShockZ2;
        vStateVariable = zeros(nStateVariable);
        vStateVariable(nIndexZ1) = log(vGridShockZ1(nShockZ1));
        vStateVariable(nIndexZ2) = log(vGridShockZ2(nShockZ2));
        for nIntangebleCapital = 1:nGridIntangebleCapital
            vStateVariable(nIndexIntangebleCapital) = vGridIntangebleCapital(nIntangebleCapital)-vSteadyStateState(nIndexIntangebleCapital);
            vPolicyIntangebleCapitalOne.IntangebleCapital(nShock,nIntangebleCapital) = vSteadyStateControl(nIndexIntangebleCapitalOne)+vFirstCoeffStateOnControl(nIndexIntangebleCapitalOne,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexIntangebleCapitalOne)...
                +mSecondCoeffStateOnControl(nIndexIntangebleCapitalOne,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexIntangebleCapitalOne,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexIntangebleCapitalOne,:)*vStateVariable;
        end
    end
end


% Changing Debt

% No shock
for nDebt = 1:nGridDebt
    vStateVariable = zeros(nStateVariable);
    vStateVariable(nIndexDebt) = vGridDebt(nDebt)-vSteadyStateState(nIndexDebt);
    vPolicyIntangebleCapitalOne.Debt(5,nDebt) = vSteadyStateControl(nIndexIntangebleCapitalOne)+vFirstCoeffStateOnControl(nIndexIntangebleCapitalOne,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexIntangebleCapitalOne)...
        +mSecondCoeffStateOnControl(nIndexIntangebleCapitalOne,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexIntangebleCapitalOne,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexIntangebleCapitalOne,:)*vStateVariable;
end

for nShockZ2 = 1:nGridShockZ1
    for nShockZ1 = 1:nGridShockZ2
        nShock = (nShockZ1-1)*nGridShockZ2+nShockZ2;
        vStateVariable = zeros(nStateVariable);
        vStateVariable(nIndexZ1) = log(vGridShockZ1(nShockZ1));
        vStateVariable(nIndexZ2) = log(vGridShockZ2(nShockZ2));
        
        for nDebt = 1:nGridDebt
            vStateVariable(nIndexDebt) = vGridDebt(nDebt)-vSteadyStateState(nIndexDebt);
            vPolicyIntangebleCapitalOne.Debt(nShock,nDebt) = vSteadyStateControl(nIndexIntangebleCapitalOne)+vFirstCoeffStateOnControl(nIndexIntangebleCapitalOne,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexIntangebleCapitalOne)...
                +mSecondCoeffStateOnControl(nIndexIntangebleCapitalOne,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexIntangebleCapitalOne,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexIntangebleCapitalOne,:)*vStateVariable;
        end
    end
end



save(resultmat,'vPolicyIntangebleCapitalOne','-append')


%%
%-------------------------------------- 
% Labor 1 Policy Function 
%---------------------------------------

% Changing Tangeble capital

% No shock
for nTangebleCapital = 1:nGridTangebleCapital
    vStateVariable = zeros(nStateVariable);
    vStateVariable(nIndexTangebleCapital) = vGridTangebleCapital(nTangebleCapital)-vSteadyStateState(nIndexTangebleCapital);
    vPolicyLaborOne.TangebleCapital(5,nTangebleCapital) = vSteadyStateControl(nIndexLaborOne)+vFirstCoeffStateOnControl(nIndexLaborOne,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexLaborOne)...
        +mSecondCoeffStateOnControl(nIndexLaborOne,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexLaborOne,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexLaborOne,:)*vStateVariable;
end

for nShockZ2 = 1:nGridShockZ1
    for nShockZ1 = 1:nGridShockZ2
        nShock = (nShockZ1-1)*nGridShockZ2+nShockZ2;
        vStateVariable = zeros(nStateVariable);
        vStateVariable(nIndexZ1) = log(vGridShockZ1(nShockZ1));
        vStateVariable(nIndexZ2) = log(vGridShockZ2(nShockZ2));
        
        for nTangebleCapital = 1:nGridTangebleCapital
            vStateVariable(nIndexTangebleCapital) = vGridTangebleCapital(nTangebleCapital)-vSteadyStateState(nIndexTangebleCapital);
            vPolicyLaborOne.TangebleCapital(nShock,nTangebleCapital) = vSteadyStateControl(nIndexLaborOne)+vFirstCoeffStateOnControl(nIndexLaborOne,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexLaborOne)...
                +mSecondCoeffStateOnControl(nIndexLaborOne,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexLaborOne,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexLaborOne,:)*vStateVariable;
        end
    end
end


% Changing Intangeble capital

% No shock
for nIntangebleCapital = 1:nGridIntangebleCapital
    vStateVariable = zeros(nStateVariable);
    vStateVariable(nIndexIntangebleCapital) = vGridIntangebleCapital(nIntangebleCapital)-vSteadyStateState(nIndexIntangebleCapital);
    vPolicyLaborOne.IntangebleCapital(5,nIntangebleCapital) = vSteadyStateControl(nIndexLaborOne)+vFirstCoeffStateOnControl(nIndexLaborOne,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexLaborOne)...
        +mSecondCoeffStateOnControl(nIndexLaborOne,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexLaborOne,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexLaborOne,:)*vStateVariable;
end

for nShockZ2 = 1:nGridShockZ1
    for nShockZ1 = 1:nGridShockZ2
        nShock = (nShockZ1-1)*nGridShockZ2+nShockZ2;
        vStateVariable = zeros(nStateVariable);
        vStateVariable(nIndexZ1) = log(vGridShockZ1(nShockZ1));
        vStateVariable(nIndexZ2) = log(vGridShockZ2(nShockZ2));
        for nIntangebleCapital = 1:nGridIntangebleCapital
            vStateVariable(nIndexIntangebleCapital) = vGridIntangebleCapital(nIntangebleCapital)-vSteadyStateState(nIndexIntangebleCapital);
            vPolicyLaborOne.IntangebleCapital(nShock,nIntangebleCapital) = vSteadyStateControl(nIndexLaborOne)+vFirstCoeffStateOnControl(nIndexLaborOne,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexLaborOne)...
                +mSecondCoeffStateOnControl(nIndexLaborOne,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexLaborOne,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexLaborOne,:)*vStateVariable;
        end
    end
end

% Changing Debt

% No shock
for nDebt = 1:nGridDebt
    vStateVariable = zeros(nStateVariable);
    vStateVariable(nIndexDebt) = vGridDebt(nDebt)-vSteadyStateState(nIndexDebt);
    vPolicyLaborOne.Debt(5,nDebt) = vSteadyStateControl(nIndexLaborOne)+vFirstCoeffStateOnControl(nIndexLaborOne,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexLaborOne)...
        +mSecondCoeffStateOnControl(nIndexLaborOne,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexLaborOne,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexLaborOne,:)*vStateVariable;
end

for nShockZ2 = 1:nGridShockZ1
    for nShockZ1 = 1:nGridShockZ2
        nShock = (nShockZ1-1)*nGridShockZ2+nShockZ2;
        vStateVariable = zeros(nStateVariable);
        vStateVariable(nIndexZ1) = log(vGridShockZ1(nShockZ1));
        vStateVariable(nIndexZ2) = log(vGridShockZ2(nShockZ2));
        
        for nDebt = 1:nGridDebt
            vStateVariable(nIndexDebt) = vGridDebt(nDebt)-vSteadyStateState(nIndexDebt);
            vPolicyLaborOne.Debt(nShock,nDebt) = vSteadyStateControl(nIndexLaborOne)+vFirstCoeffStateOnControl(nIndexLaborOne,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexLaborOne)...
                +mSecondCoeffStateOnControl(nIndexLaborOne,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexLaborOne,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexLaborOne,:)*vStateVariable;
        end
    end
end


save(resultmat,'vPolicyLaborOne','-append')


%%
%-------------------------------------- 
% Labor 2 Policy Function 
%---------------------------------------

% Changing Tangeble capital

% No shock
for nTangebleCapital = 1:nGridTangebleCapital
    vStateVariable = zeros(nStateVariable);
    vStateVariable(nIndexTangebleCapital) = vGridTangebleCapital(nTangebleCapital)-vSteadyStateState(nIndexTangebleCapital);
    vPolicyLaborTwo.TangebleCapital(5,nTangebleCapital) = vSteadyStateControl(nIndexLaborTwo)+vFirstCoeffStateOnControl(nIndexLaborTwo,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexLaborTwo)...
        +mSecondCoeffStateOnControl(nIndexLaborTwo,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexLaborTwo,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexLaborTwo,:)*vStateVariable;
end

for nShockZ2 = 1:nGridShockZ1
    for nShockZ1 = 1:nGridShockZ2
        nShock = (nShockZ1-1)*nGridShockZ2+nShockZ2;
        vStateVariable = zeros(nStateVariable);
        vStateVariable(nIndexZ1) = log(vGridShockZ1(nShockZ1));
        vStateVariable(nIndexZ2) = log(vGridShockZ2(nShockZ2));
        
        for nTangebleCapital = 1:nGridTangebleCapital
            vStateVariable(nIndexTangebleCapital) = vGridTangebleCapital(nTangebleCapital)-vSteadyStateState(nIndexTangebleCapital);
            vPolicyLaborTwo.TangebleCapital(nShock,nTangebleCapital) = vSteadyStateControl(nIndexLaborTwo)+vFirstCoeffStateOnControl(nIndexLaborTwo,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexLaborTwo)...
                +mSecondCoeffStateOnControl(nIndexLaborTwo,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexLaborTwo,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexLaborTwo,:)*vStateVariable;
        end
    end
end


% Changing Intangeble capital

% No shock
for nIntangebleCapital = 1:nGridIntangebleCapital
    vStateVariable = zeros(nStateVariable);
    vStateVariable(nIndexIntangebleCapital) = vGridIntangebleCapital(nIntangebleCapital)-vSteadyStateState(nIndexIntangebleCapital);
    vPolicyLaborTwo.IntangebleCapital(5,nIntangebleCapital) = vSteadyStateControl(nIndexLaborTwo)+vFirstCoeffStateOnControl(nIndexLaborTwo,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexLaborTwo)...
        +mSecondCoeffStateOnControl(nIndexLaborTwo,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexLaborTwo,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexLaborTwo,:)*vStateVariable;
end

for nShockZ2 = 1:nGridShockZ1
    for nShockZ1 = 1:nGridShockZ2
        nShock = (nShockZ1-1)*nGridShockZ2+nShockZ2;
        vStateVariable = zeros(nStateVariable);
        vStateVariable(nIndexZ1) = log(vGridShockZ1(nShockZ1));
        vStateVariable(nIndexZ2) = log(vGridShockZ2(nShockZ2));
        for nIntangebleCapital = 1:nGridIntangebleCapital
            vStateVariable(nIndexIntangebleCapital) = vGridIntangebleCapital(nIntangebleCapital)-vSteadyStateState(nIndexIntangebleCapital);
            vPolicyLaborTwo.IntangebleCapital(nShock,nIntangebleCapital) = vSteadyStateControl(nIndexLaborTwo)+vFirstCoeffStateOnControl(nIndexLaborTwo,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexLaborTwo)...
                +mSecondCoeffStateOnControl(nIndexLaborTwo,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexLaborTwo,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexLaborTwo,:)*vStateVariable;
        end
    end
end

% Changing Debt

% No shock
for nDebt = 1:nGridDebt
    vStateVariable = zeros(nStateVariable);
    vStateVariable(nIndexDebt) = vGridDebt(nDebt)-vSteadyStateState(nIndexDebt);
    vPolicyLaborTwo.Debt(5,nDebt) = vSteadyStateControl(nIndexLaborTwo)+vFirstCoeffStateOnControl(nIndexLaborTwo,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexLaborTwo)...
        +mSecondCoeffStateOnControl(nIndexLaborTwo,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexLaborTwo,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexLaborTwo,:)*vStateVariable;
end

for nShockZ2 = 1:nGridShockZ1
    for nShockZ1 = 1:nGridShockZ2
        nShock = (nShockZ1-1)*nGridShockZ2+nShockZ2;
        vStateVariable = zeros(nStateVariable);
        vStateVariable(nIndexZ1) = log(vGridShockZ1(nShockZ1));
        vStateVariable(nIndexZ2) = log(vGridShockZ2(nShockZ2));
        
        for nDebt = 1:nGridDebt
            vStateVariable(nIndexDebt) = vGridDebt(nDebt)-vSteadyStateState(nIndexDebt);
            vPolicyLaborTwo.Debt(nShock,nDebt) = vSteadyStateControl(nIndexLaborTwo)+vFirstCoeffStateOnControl(nIndexLaborTwo,:)*vStateVariable+mSecondCoeffPerturbationOnControl(nIndexLaborTwo)...
                +mSecondCoeffStateOnControl(nIndexLaborTwo,:)*kron(vStateVariable,vStateVariable)+mThridCoeffStateOnControl(nIndexLaborTwo,:)*kron(kron(vStateVariable,vStateVariable),vStateVariable)+mThirdCoeffCrossTermOnControl(nIndexLaborTwo,:)*vStateVariable;
        end
    end
end


save(resultmat,'vPolicyLaborTwo','-append')

%% Plot 


figure1 = figure(1);
name1 = strcat('ValueFunction_Perturbation.pdf');

subplot(2,2,1)
plot(vGridTangebleCapital',vValueFunction.TangebleCapital','LineWidth',lwidth)
axis tight;
legend('both high','z1 high','z2 high','both low','No shock')
set(gca,'FontSize',fsize)
title('Value Function')
xlabel('Tangeble Capital','FontSize',fsize)
set(legend,'FontSize',fsize1)
set(legend,'Position',[0.22,0.85,0.01,0.01])
set(legend,'Box','off')


subplot(2,2,2)
plot(vGridIntangebleCapital',vValueFunction.IntangebleCapital','LineWidth',lwidth)
axis tight;
set(gca,'FontSize',fsize)
title('Value Function')
xlabel('Intangeble Capital','FontSize',fsize)


subplot(2,2,3)
plot(vGridDebt',vValueFunction.Debt','LineWidth',lwidth)
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
name2 = strcat('PolicyFunction_TangebleCapital_Perturbation.pdf');

subplot(2,2,1)
plot(vGridTangebleCapital',vPolicyTangebleCapital.TangebleCapital','LineWidth',lwidth)
axis tight;
legend('both high','z1 high','z2 high','both low','No shock')
title('Tangeble Capital Pol Fun')
xlabel('Tangeble Capital','FontSize',fsize)
set(legend,'FontSize',fsize1)
set(legend,'Position',[0.22,0.85,0.01,0.01])
set(legend,'Box','off')
set(gca,'FontSize',fsize)

subplot(2,2,2)
plot(vGridIntangebleCapital',vPolicyTangebleCapital.IntangebleCapital','LineWidth',lwidth)
axis tight;
title('Tangeble Capital Pol Fun')
xlabel('Intangeble Capital','FontSize',fsize)
set(gca,'FontSize',fsize)

subplot(2,2,3)
plot(vGridDebt',vPolicyTangebleCapital.Debt','LineWidth',lwidth)
axis tight;
title('Tangeble Capital Pol Fun')
xlabel('Debt','FontSize',fsize)
set(gca,'FontSize',fsize)


cd(path_Figure);
set(gcf,'PaperOrientation','landscape','PaperPosition',[0 0 16 12])
set(gcf, 'PaperSize', [16 12]);
print(figure2,'-dpdf',name2)
cd(path_Main);


figure3 = figure(3);
name3 = strcat('PolicyFunction_IntangebleCapital_Perturbation.pdf');
subplot(2,2,1)
plot(vGridTangebleCapital',vPolicyIntangebleCapital.TangebleCapital','LineWidth',lwidth)
axis tight;
legend('both high','z1 high','z2 high','both low','No shock')
title('Intangeble Capital Pol Fun')
xlabel('Tangeble Capital','FontSize',fsize)
set(legend,'FontSize',fsize1)
set(legend,'Position',[0.22,0.85,0.01,0.01])
set(legend,'Box','off')
set(gca,'FontSize',fsize)



subplot(2,2,2)
plot(vGridIntangebleCapital',vPolicyIntangebleCapital.IntangebleCapital','LineWidth',lwidth)
axis tight;
title('Intangeble Capital Pol Fun')
xlabel('Intangeble Capital','FontSize',fsize)
set(gca,'FontSize',fsize)


subplot(2,2,3)
plot(vGridDebt',vPolicyIntangebleCapital.Debt','LineWidth',lwidth)
axis tight;
title('Intangeble Capital Pol Fun')
xlabel('Debt','FontSize',fsize)
set(gca,'FontSize',fsize)


cd(path_Figure);
set(gcf,'PaperOrientation','landscape','PaperPosition',[0 0 16 12])
set(gcf, 'PaperSize', [16 12]);
print(figure3,'-dpdf',name3)
cd(path_Main);



figure4 = figure(4);
name4=strcat('PolicyFunction_TangebleCapitalOne_Perturbation.pdf');
subplot(2,2,1)
plot(vGridTangebleCapital',vPolicyTangebleCapitalOne.TangebleCapital','LineWidth',lwidth)
axis tight;
legend('both high','z1 high','z2 high','both low','No shock')
title('Tangeble Capital 1 Pol Fun')
xlabel('Tangeble Capital','FontSize',fsize)
set(legend,'FontSize',fsize1)
set(legend,'Position',[0.22,0.85,0.01,0.01])
set(legend,'Box','off')
set(gca,'FontSize',fsize)



subplot(2,2,2)
plot(vGridIntangebleCapital',vPolicyTangebleCapitalOne.IntangebleCapital','LineWidth',lwidth)
axis tight;
title('Tangeble Capital 1 Pol Fun')
xlabel('Intangeble Capital','FontSize',fsize)
set(gca,'FontSize',fsize)

subplot(2,2,3)
plot(vGridDebt',vPolicyTangebleCapitalOne.Debt','LineWidth',lwidth)
axis tight;
title('Tangeble Capital 1 Pol Fun')
xlabel('Debt','FontSize',fsize)
set(gca,'FontSize',fsize)



cd(path_Figure);
set(gcf,'PaperOrientation','landscape','PaperPosition',[0 0 16 12])
set(gcf, 'PaperSize', [16 12]);
print(figure4,'-dpdf',name4)
cd(path_Main);



figure5 = figure(5);
name5=strcat('PolicyFunction_IntangebleCapitalOne_Perturbation.pdf');
subplot(2,2,1)
plot(vGridTangebleCapital',vPolicyIntangebleCapitalOne.TangebleCapital','LineWidth',lwidth)
axis tight;
title('Intangeble Capital 1 Pol Fun')
legend('both high','z1 high','z2 high','both low','No shock')
xlabel('Tangeble Capital','FontSize',fsize)
set(legend,'FontSize',fsize1)
set(legend,'Position',[0.22,0.82,0.01,0.01])
set(legend,'Box','off')
set(gca,'FontSize',fsize)



subplot(2,2,2)
plot(vGridIntangebleCapital',vPolicyIntangebleCapitalOne.IntangebleCapital','LineWidth',lwidth)
axis tight;
title('Intangeble Capital 1 Pol Fun')
xlabel('Intangeble Capital','FontSize',fsize)
set(gca,'FontSize',fsize)



subplot(2,2,3)
plot(vGridDebt',vPolicyIntangebleCapitalOne.Debt','LineWidth',lwidth)
axis tight;
title('Intangeble Capital 1 Pol Fun')
xlabel('Debt','FontSize',fsize)
set(gca,'FontSize',fsize)

cd(path_Figure);
set(gcf,'PaperOrientation','landscape','PaperPosition',[0 0 16 12])
set(gcf, 'PaperSize', [16 12]);
print(figure5,'-dpdf',name5)
cd(path_Main);




figure6 = figure(6);
name6=strcat('PolicyFunction_LaborOne_Perturbation.pdf');
subplot(2,2,1)
plot(vGridTangebleCapital',vPolicyLaborOne.TangebleCapital','LineWidth',lwidth)
axis tight;
legend('both high','z1 high','z2 high','both low','No shock')
title('Labor 1 Pol Fun')
xlabel('Tangeble Capital','FontSize',fsize)
set(legend,'FontSize',fsize1)
set(legend,'Position',[0.22,0.85,0.01,0.01])
set(legend,'Box','off')
set(gca,'FontSize',fsize)



subplot(2,2,2)
plot(vGridIntangebleCapital',vPolicyLaborOne.IntangebleCapital','LineWidth',lwidth)
axis tight;
title('Labor 1 Pol Fun')
xlabel('Intangeble Capital','FontSize',fsize)
set(gca,'FontSize',fsize)

subplot(2,2,3)
plot(vGridDebt',vPolicyLaborOne.Debt','LineWidth',lwidth)
axis tight;
title('Labor 1 Pol Fun')
xlabel('Debt','FontSize',fsize)
set(gca,'FontSize',fsize)


cd(path_Figure);
set(gcf,'PaperOrientation','landscape','PaperPosition',[0 0 16 12])
set(gcf, 'PaperSize', [16 12]);
print(figure6,'-dpdf',name6)
cd(path_Main);

figure7 = figure(7);
name7=strcat('PolicyFunction_LaborTwo_Perturbation.pdf');
subplot(2,2,1)
plot(vGridTangebleCapital',vPolicyLaborTwo.TangebleCapital','LineWidth',lwidth)
axis tight;
legend('both high','z1 high','z2 high','both low','No shock')
title('Labor 2 Pol Fun')
xlabel('Tangeble Capital','FontSize',fsize)
set(legend,'FontSize',fsize1)
set(legend,'Position',[0.22,0.85,0.01,0.01])
set(legend,'Box','off')
set(gca,'FontSize',fsize)



subplot(2,2,2)
plot(vGridIntangebleCapital',vPolicyLaborTwo.IntangebleCapital','LineWidth',lwidth)
axis tight;
title('Labor 2 Pol Fun')
xlabel('Intangeble Capital','FontSize',fsize)
set(gca,'FontSize',fsize)

subplot(2,2,3)
plot(vGridDebt',vPolicyLaborTwo.Debt','LineWidth',lwidth)
axis tight;
title('Labor 2 Pol Fun')
xlabel('Debt','FontSize',fsize)
set(gca,'FontSize',fsize)

cd(path_Figure);
set(gcf,'PaperOrientation','landscape','PaperPosition',[0 0 16 12])
set(gcf, 'PaperSize', [16 12]);
print(figure7,'-dpdf',name7)
cd(path_Main);


figure8 = figure(8);
name8=strcat('PolicyFunction_Debt_Perturbation.pdf');
subplot(2,2,1)
plot(vGridTangebleCapital',vPolicyDebt.TangebleCapital','LineWidth',lwidth)
axis tight;
legend('both high','z1 high','z2 high','both low','No shock')
title('Debt Pol Fun')
xlabel('Tangeble Capital','FontSize',fsize)
set(legend,'FontSize',fsize1)
set(legend,'Position',[0.22,0.85,0.01,0.01])
set(legend,'Box','off')
set(gca,'FontSize',fsize)



subplot(2,2,2)
plot(vGridIntangebleCapital',vPolicyDebt.IntangebleCapital','LineWidth',lwidth)
axis tight;
title('Debt Pol Fun')
xlabel('Intangeble Capital','FontSize',fsize)
set(gca,'FontSize',fsize)

subplot(2,2,3)
plot(vGridDebt',vPolicyDebt.Debt','LineWidth',lwidth)
axis tight;
title('Debt Pol Fun')
xlabel('Debt','FontSize',fsize)
set(gca,'FontSize',fsize)

cd(path_Figure);
set(gcf,'PaperOrientation','landscape','PaperPosition',[0 0 16 12])
set(gcf, 'PaperSize', [16 12]);
print(figure8,'-dpdf',name8)
cd(path_Main);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  3. Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nGridPeriod = optPruning.numSim;
vGridPeriod = 1:1:nGridPeriod;
% Because we already have the simulated data, we simply plot them

name99 = strcat('Dynamics_state_Perturbation.pdf');
figure99 = figure(99);
subplot(3,2,1)
vPathTangebleCapial = outDynare.vSim(nIndexTangebleCapital,:);
plot(vGridPeriod,vPathTangebleCapial)
title('Tangeble Capital')
axis tight;
set(gca,'FontSize',fsize)

subplot(3,2,2)
vPathIntangebleCapial = outDynare.vSim(nIndexIntangebleCapital,:);
plot(vGridPeriod,vPathIntangebleCapial)
title('Intangeble Capital')
axis tight;
set(gca,'FontSize',fsize)

subplot(3,2,3)
vPathDebt = outDynare.vSim(nIndexDebt,:);
plot(vGridPeriod,vPathDebt)
title('Debt')
axis tight;
set(gca,'FontSize',fsize)

subplot(3,2,4)
vPathez1 = outDynare.vSim(nIndexZ1,:);
plot(vGridPeriod,vPathez1)
title('z1')
axis tight;
set(gca,'FontSize',fsize)

subplot(3,2,5)
vPathez2 = outDynare.vSim(nIndexZ2,:);
plot(vGridPeriod,vPathez2)
title('z2')
axis tight;
set(gca,'FontSize',fsize)

cd(path_Figure);
set(gcf,'PaperOrientation','landscape','PaperPosition',[0 0 15 10])
set(gcf, 'PaperSize', [15 10]);
print(figure99,'-dpdf',name99)
cd(path_Main);


name100 = strcat('Dynamics_endo_Perturbation.pdf');
figure100 = figure(100);
subplot(3,2,1)
vPathLaborOne = outDynare.ySim(nIndexLaborOne,:);
plot(vGridPeriod,vPathLaborOne)
title('Labor 1')
axis tight;
set(gca,'FontSize',fsize)

subplot(3,2,2)
vPathLaborTwo = outDynare.ySim(nIndexLaborTwo,:);
plot(vGridPeriod,vPathLaborTwo)
title('Labor 2')
axis tight;
set(gca,'FontSize',fsize)

subplot(3,2,3)
vPathTangebleCapitalOne = outDynare.ySim(nIndexTangebleCapitalOne,:);
plot(vGridPeriod,vPathTangebleCapitalOne)
title('Tangeble Capital 1')
axis tight;
set(gca,'FontSize',fsize)

subplot(3,2,4)
vPathIntangebleCapitalOne = outDynare.ySim(nIndexIntangebleCapitalOne,:);
plot(vGridPeriod,vPathIntangebleCapitalOne)
title('Intangeble Capital 1')
axis tight;
set(gca,'FontSize',fsize)

subplot(3,2,5)
vPathValue = outDynare.ySim(nIndexValueFunction,:);
plot(vGridPeriod,vPathValue)
title('Value Function')
axis tight;
set(gca,'FontSize',fsize)





cd(path_Figure);
set(gcf,'PaperOrientation','landscape','PaperPosition',[0 0 15 10])
set(gcf, 'PaperSize', [15 10]);
print(figure100,'-dpdf',name100)
cd(path_Main);


