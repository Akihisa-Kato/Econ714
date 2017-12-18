function [mTheta] = Residual_Function(mTheta)

% Residual_Function
%   Compute the residuals of the optimality conditions in Econ714 HW2 Q1
%   This function will be used inside of fsolve to find Smolyak-Chebychev
%   polynomial coefficients mTheta
%
%   INPUT: mTheta : nVariable*G(5,3)*nGridProductivityI*nGridProductivityI matrix 
%          Variable is ordered as 
%          [J kt kt1 ki1 l1 l2 b mmu2]';
%
%   OUTPUT: Diff : differences bewteen RHS and LHS of optimality conditions
%   nVariable*G(5,3)*nGridProductivityI*nGridProductivityI matrix

global llambda1 llambda2 ttheta1 ttheta2 aalpha1 aalpha2 ggamma1 ggamma2 ddeltat ddeltai w ri r nVariable ...
    vProductivityT vProductivityI mTransitionT mTransitionI nMaxOrder nGridProductivityI nGridProductivityT ...
    nIndexJ nIndexkt nIndexl1 nIndexl2 nIndexki2 nIndexkt2 nIndexb nIndexki...
    nTangebleCapitalMax nTangebleCapitalMin nIntangebleCapitalMax nIntangebleCapitalMin nDebtMax nDebtMin mSmolyakGrid...
    nGridSmolyak vChebychev_kt vChebychev_ki vChebychev_b nGridProductivity


mmu1 = 1;

options=optimset('Display','Iter','TolFun',10^(-15),'TolX',10^(-15),'FunValCheck','on','MaxFunEval',800^100,'MaxIter',10^10);%,'MaxIter',6);
mTheta = fsolve(@notime_iter,mTheta,options);

function Diff = notime_iter(mTheta)
%Diff_sec = zeros(nGridSmolyak*nVariable,1);

Diff = zeros(nVariable,nGridSmolyak,nGridProductivityI*nGridProductivityI);

flag = 0;
for nProductivityT = 1:nGridProductivityT

    ZT = vProductivityT(nProductivityT);

    for nProductivityI = 1:nGridProductivityI

        ZI = vProductivityI(nProductivityI);
        nProductivity = (nProductivityT-1)*nGridProductivityT+nProductivityI;
        
        vTangebleCapitalTomorrow = zeros(nGridSmolyak,1);
        vIntangebleCapitalTomorrow = zeros(nGridSmolyak,1);
        vDebtTomorrow = zeros(nGridSmolyak,1);
        vTangebleCapitalTwo = zeros(nGridSmolyak,1);
        vIntangebleCapitalTwo = zeros(nGridSmolyak,1);
        vLaborOne = zeros(nGridSmolyak,1);
        vLaborTwo = zeros(nGridSmolyak,1);
        vValue = zeros(nGridSmolyak,1);
        vMmuTwo = zeros(nGridSmolyak,1);
        vg1 = zeros(nGridSmolyak,1);
        vg2 = zeros(nGridSmolyak,1);
        vInvestment = zeros(nGridSmolyak,1);
        vSale = zeros(nGridSmolyak,1);
        vDividend = zeros(nGridSmolyak,1);

        for nStateGrid = 1:nGridSmolyak % Sparse grids
            
            vChebychevToday = [1; vChebychev_kt(nStateGrid,2:nMaxOrder)'; vChebychev_ki(nStateGrid,2:nMaxOrder)'; vChebychev_b(nStateGrid,2:nMaxOrder)'; % const and no cross terms (13)
                               vChebychev_kt(nStateGrid,2)*vChebychev_ki(nStateGrid,2); vChebychev_kt(nStateGrid,2)*vChebychev_b(nStateGrid,2); vChebychev_ki(nStateGrid,2)*vChebychev_b(nStateGrid,2); % cross 1*1 (3)
                               vChebychev_kt(nStateGrid,2)*vChebychev_ki(nStateGrid,3); vChebychev_kt(nStateGrid,2)*vChebychev_b(nStateGrid,3); vChebychev_ki(nStateGrid,2)*vChebychev_kt(nStateGrid,3); 
                               vChebychev_ki(nStateGrid,2)*vChebychev_b(nStateGrid,3);  vChebychev_b(nStateGrid,2)*vChebychev_kt(nStateGrid,3); vChebychev_b(nStateGrid,2)*vChebychev_ki(nStateGrid,3); % cross 1*2 (6)
                               vChebychev_kt(nStateGrid,3)*vChebychev_ki(nStateGrid,3); vChebychev_kt(nStateGrid,3)*vChebychev_b(nStateGrid,3); vChebychev_ki(nStateGrid,3)*vChebychev_b(nStateGrid,3)];% cross 2*2 (3)


            nTangebleCapitalToday = (mSmolyakGrid(nStateGrid,1)+1)*(nTangebleCapitalMax-nTangebleCapitalMin)/2+nTangebleCapitalMin;
            nIntangebleCapitalToday = (mSmolyakGrid(nStateGrid,2)+1)*(nIntangebleCapitalMax-nIntangebleCapitalMin)/2+nIntangebleCapitalMin;
            nDebtToday = (mSmolyakGrid(nStateGrid,3)+1)*(nDebtMax-nDebtMin)/2+nDebtMin;

            


                    % Policy function for kt1 ki1 l1 l2 kt' b' ki'
                    % 1. kt'
                    vTangebleCapitalTomorrow(nStateGrid) = vec(mTheta(nIndexkt,:,nProductivity))'*vChebychevToday;
                    %{
                    if (vTangebleCapitalTomorrow(nStateGrid) < nTangebleCapitalMin)
                        vTangebleCapitalTomorrow(nStateGrid) = nTangebleCapitalMin+0.01;
                        %disp('TangebleCapitalTomorrow break lower bound')
                        %flag = 1;
                    elseif (vTangebleCapitalTomorrow(nStateGrid) > nTangebleCapitalMax)
                        vTangebleCapitalTomorrow(nStateGrid) = nTangebleCapitalMax-0.01;
                        %disp('TangebleCapitalTomorrow break upper bound')
                        %flag = 1;
                    end
                    %}

                    % 2. kt1 
                    vTangebleCapitalTwo(nStateGrid) = vec(mTheta(nIndexkt2,:,nProductivity))'*vChebychevToday;

                    % 3. kti
                    vIntangebleCapitalTwo(nStateGrid) = vec(mTheta(nIndexki2,:,nProductivity))'*vChebychevToday;

                    % 4. l1
                    vLaborOne(nStateGrid) = vec(mTheta(nIndexl1,:,nProductivity))'*vChebychevToday;

                    % 5. l2
                    vLaborTwo(nStateGrid) = vec(mTheta(nIndexl2,:,nProductivity))'*vChebychevToday;

                    % 6. b'
                    vDebtTomorrow(nStateGrid) = vec(mTheta(nIndexb,:,nProductivity))'*vChebychevToday;
                    %{
                    if (vDebtTomorrow(nStateGrid) < nDebtMin)
                        vDebtTomorrow(nStateGrid) = nDebtMin+0.01;
                        %disp('DebtTomorrow break lower bound')
                        %flag = 3;
                    elseif (vDebtTomorrow(nStateGrid) > nDebtMax)
                        vDebtTomorrow(nStateGrid) = nDebtMax-0.01;
                        %disp('DebtTomorrow break upper bound')
                        %flag = 3;
                    end
                    %}

                    % 7. mmu2
                    %vMmuTwo(nStateGrid) = vec(mTheta(nIndexmmu2,:,nProductivity))'*vChebychevToday;
                    
                    
                    
                    % 8 & 9 Construct auxiliary variables (g1 g2)

                    vg1(nStateGrid) = ttheta1*(nTangebleCapitalToday-vTangebleCapitalTwo(nStateGrid))^((llambda1-1)/llambda1)+(1-ttheta1)*(nIntangebleCapitalToday-vIntangebleCapitalTwo(nStateGrid))^((llambda1-1)/llambda1);
                    vg2(nStateGrid) = ttheta2*vTangebleCapitalTwo(nStateGrid)^((llambda2-1)/llambda2)+(1-ttheta2)*vIntangebleCapitalTwo(nStateGrid)^((llambda2-1)/llambda2);

                    % 9. ki'
                    %vIntangebleCapitalTomorrow(nStateGrid) = vg2(nStateGrid)^(aalpha2*llambda2/(llambda2-1))*(exp(ZI)*vLaborTwo(nStateGrid))^ggamma2+(1-ddeltai)*nIntangebleCapitalToday;
                    vIntangebleCapitalTomorrow(nStateGrid) = vec(mTheta(nIndexki,:,nProductivity))'*vChebychevToday;
                    %{
                    if (vIntangebleCapitalTomorrow(nStateGrid) < nIntangebleCapitalMin)
                        vIntangebleCapitalTomorrow(nStateGrid) = nIntangebleCapitalMin+0.01;
                        %disp('IntangebleCapitalTomorrow break lower bound')
                        %flag = 2;
                    elseif (vIntangebleCapitalTomorrow(nStateGrid) > nIntangebleCapitalMax)
                        vIntangebleCapitalTomorrow(nStateGrid) = nIntangebleCapitalMax-0.01;
                        %disp('IntangebleCapitalTomorrow break upper bound')
                        %flag = 2;
                    end
                    %}
                    % 10. x, investment
                    
                    vInvestment(nStateGrid) = vTangebleCapitalTomorrow(nStateGrid)-(1-ddeltat)*nTangebleCapitalToday;
                    
                    % 11. s, sales
                    vSale(nStateGrid) = vg1(nStateGrid)^(aalpha1*llambda1/(llambda1-1))*(exp(ZT)*vLaborOne(nStateGrid))^ggamma1-vInvestment(nStateGrid);
                    
                    % 12. d, Dividend
                    vDividend(nStateGrid) = vSale(nStateGrid)-w*(vLaborOne(nStateGrid)+vLaborTwo(nStateGrid))+vDebtTomorrow(nStateGrid)-(1+r)*nDebtToday...
                        -0.02*(vDebtTomorrow(nStateGrid)-nDebtToday)^2-0.01*(nDebtToday-0.2)^2;
                    
                    
                    % ki1
                    vTangebleCapitalOne(nStateGrid,1) = nTangebleCapitalToday-vTangebleCapitalTwo(nStateGrid);
                    vIntangebleCapitalOne(nStateGrid,1) = nIntangebleCapitalToday-vIntangebleCapitalTwo(nStateGrid);
                    vDebtChange(nStateGrid,1) = vDebtTomorrow(nStateGrid)-nDebtToday;
                    



                    % Value Function
                    vValue(nStateGrid) = vec(mTheta(nIndexJ,:,nProductivity))'*vChebychevToday;
                    
                    
                    

        end
    
    
        % Scale Tomorrow's state into [-1,1]
        % nTangebleCapitalTomorrow
        vTangebleCapitalTomorrowScaleDown = 2*((vTangebleCapitalTomorrow-nTangebleCapitalMin)/(nTangebleCapitalMax-nTangebleCapitalMin))-1;

        % nIntangebleCapitalTomorrow
        vIntangebleCapitalTomorrowScaleDown = 2*((vIntangebleCapitalTomorrow-nIntangebleCapitalMin)/(nIntangebleCapitalMax-nIntangebleCapitalMin))-1;

        % nDebtTomorrow
        vDebtTomorrowScaleDown = 2*((vDebtTomorrow-nDebtMin)/(nDebtMax-nDebtMin))-1;


        % value of polynomials at each scaled tomorrow's state (scaled down)
        vChebychev_kt_prime = ones(nGridSmolyak,nMaxOrder);
        vChebychev_ki_prime = ones(nGridSmolyak,nMaxOrder);
        vChebychev_b_prime  = ones(nGridSmolyak,nMaxOrder);
        vChebychev_kt_prime(:,2) = vTangebleCapitalTomorrowScaleDown;
        vChebychev_ki_prime(:,2) = vIntangebleCapitalTomorrowScaleDown;
        vChebychev_b_prime(:,2) = vDebtTomorrowScaleDown;

        % compute values of chebychev pols
        for nOrder=3:nMaxOrder
            vChebychev_kt_prime(:,nOrder) = 2*vChebychev_kt_prime(:,2).*vChebychev_kt_prime(:,nOrder-1)-vChebychev_kt_prime(:,nOrder-2);
            vChebychev_ki_prime(:,nOrder) = 2*vChebychev_ki_prime(:,2).*vChebychev_ki_prime(:,nOrder-1)-vChebychev_ki_prime(:,nOrder-2);
            vChebychev_b_prime(:,nOrder)  = 2*vChebychev_b_prime(:,2).*vChebychev_b_prime(:,nOrder-1)-vChebychev_b_prime(:,nOrder-2);
        end

        
                          
        % NExt loop will be computing Expected terms
        vValuep = zeros(nGridProductivity,1);
        vJktp = zeros(nGridProductivity,1);
        vJbp = zeros(nGridProductivity,1);
        vJkip = zeros(nGridProductivity,1);
        
        
        
        % Calculate residual
        for nStateGrid = 1:nGridSmolyak % Sparse grids
            vChebychevTomorrow = [1; vChebychev_kt_prime(nStateGrid,2:nMaxOrder)'; vChebychev_ki_prime(nStateGrid,2:nMaxOrder)'; vChebychev_b_prime(nStateGrid,2:nMaxOrder)'; % const and no cross terms (13)
                                  vChebychev_kt_prime(nStateGrid,2)*vChebychev_ki_prime(nStateGrid,2); vChebychev_kt_prime(nStateGrid,2)*vChebychev_b_prime(nStateGrid,2); vChebychev_ki_prime(nStateGrid,2)*vChebychev_b_prime(nStateGrid,2); % cross 1*1 (3)
                                  vChebychev_kt_prime(nStateGrid,2)*vChebychev_ki_prime(nStateGrid,3); vChebychev_kt_prime(nStateGrid,2)*vChebychev_b_prime(nStateGrid,3); vChebychev_ki_prime(nStateGrid,2)*vChebychev_kt_prime(nStateGrid,3); 
                                  vChebychev_ki_prime(nStateGrid,2)*vChebychev_b_prime(nStateGrid,3);  vChebychev_b_prime(nStateGrid,2)*vChebychev_kt_prime(nStateGrid,3); vChebychev_b_prime(nStateGrid,2)*vChebychev_ki_prime(nStateGrid,3); % cross 1*2 (6)
                                  vChebychev_kt_prime(nStateGrid,3)*vChebychev_ki_prime(nStateGrid,3); vChebychev_kt_prime(nStateGrid,3)*vChebychev_b_prime(nStateGrid,3); vChebychev_ki_prime(nStateGrid,3)*vChebychev_b_prime(nStateGrid,3)];% cross 2*2 (3)

            %nTangebleCapitalToday = (mSmolyakGrid(nStateGrid,1)+1)*(nTangebleCapitalMax-nTangebleCapitalMin)/2+nTangebleCapitalMin;
            %nIntangebleCapitalToday = (mSmolyakGrid(nStateGrid,2)+1)*(nIntangebleCapitalMax-nIntangebleCapitalMin)/2+nIntangebleCapitalMin;
            %nDebtToday = (mSmolyakGrid(nStateGrid,3)+1)*(nDebtMax-nDebtMin)/2+nDebtMin;
            
            for nProductivityTp = 1:nGridProductivityT

                ZTp = vProductivityT(nProductivityTp);

                for nProductivityIp = 1:nGridProductivityI

                    ZIp = vProductivityI(nProductivityIp);
                    nProductivityp = (nProductivityTp-1)*nGridProductivityT+nProductivityIp;

                    % 1.Value Function Tomorrow
                    vValuep(nProductivityp) = vec(mTheta(nIndexJ,:,nProductivityp))'*vChebychevTomorrow;

                    % 2. kt1 Tomorrow
                    nTangebleCapitalTwop = vec(mTheta(nIndexkt2,:,nProductivityp))'*vChebychevTomorrow;

                    % 3. kti Tomorrow
                    nIntangebleCapitalTwop = vec(mTheta(nIndexki2,:,nProductivityp))'*vChebychevTomorrow;

                    % 4. l1 Tomorrow
                    nLaborOnep = vec(mTheta(nIndexl1,:,nProductivityp))'*vChebychevTomorrow;

                    % 5. l2 Tomorrow
                    nLaborTwop = vec(mTheta(nIndexl2,:,nProductivityp))'*vChebychevTomorrow;

                    % 6. b' Tomorrow
                    nDebtTomorrowp = vec(mTheta(nIndexb,:,nProductivityp))'*vChebychevTomorrow;
                    
                    
                    % g1 and g2
                    ng1p = ttheta1*(vTangebleCapitalTomorrow(nStateGrid)-nTangebleCapitalTwop)^((llambda1-1)/llambda1)+(1-ttheta1)*(vIntangebleCapitalTomorrow(nStateGrid)-nIntangebleCapitalTwop)^((llambda1-1)/llambda1);
                    ng2p = ttheta2*nTangebleCapitalTwop^((llambda2-1)/llambda2)+(1-ttheta2)*nIntangebleCapitalTwop^((llambda2-1)/llambda2);

                    
                    % 7. mmu2
                    %nMmu2p = vec(mTheta(nIndexmmu2,:,nProductivityp))'*vChebychevTomorrow;
                    nMmu2p = w*nLaborTwop*(ng2p^(aalpha2*llambda2/(llambda2-1))*ggamma2*(exp(ZIp)*nLaborTwop)^ggamma2)^-1;
                    % 8. Jb tomorrow
                    vJbp(nProductivityp) = -(1+r)+0.04*(nDebtTomorrowp-vDebtTomorrow(nStateGrid))-0.02*(vDebtTomorrow(nStateGrid)-0.2);
                    
                    
                    % 9. Jkt tomorrow
                    vJktp(nProductivityp) = mmu1*(1-ddeltat)+mmu1*aalpha1*ng1p^(llambda1/(llambda1-1)*aalpha1-1)*(exp(ZTp)*nLaborOnep)^ggamma1*ttheta1*(vTangebleCapitalTomorrow(nStateGrid)-nTangebleCapitalTwop)^(-1/llambda1);
                    
                    % 10. Jki tomorrow
                    vJkip(nProductivityp) = nMmu2p*(1-ddeltai)+mmu1*aalpha1*ng1p^(llambda1/(llambda1-1)*aalpha1-1)*(exp(ZTp)*nLaborOnep)^ggamma1*(1-ttheta1)*(vIntangebleCapitalTomorrow(nStateGrid)-nIntangebleCapitalTwop)^(-1/llambda1);
                    
                    
                end
            end
            
            % compute expected term    
            EJ = kron(mTransitionT(nProductivityT,:),mTransitionI(nProductivityI,:))*vValuep;
            EJkt = kron(mTransitionT(nProductivityT,:),mTransitionI(nProductivityI,:))*vJktp;
            EJki = kron(mTransitionT(nProductivityT,:),mTransitionI(nProductivityI,:))*vJkip;
            EJb = kron(mTransitionT(nProductivityT,:),mTransitionI(nProductivityI,:))*vJbp;

            
            vMmuTwo(nStateGrid) = EJki/(1+ri);
            
            
            % Residuals
            Diff(nIndexJ,nStateGrid,nProductivity) = -vValue(nStateGrid)+vDividend(nStateGrid)+EJ/(1+ri);
            
            Diff(nIndexkt,nStateGrid,nProductivity) = mmu1-EJkt/(1+ri);
            
            LHSkt2 = mmu1*aalpha1*vg1(nStateGrid)^(aalpha1*llambda1/(llambda1-1)-1)*(exp(ZT)*vLaborOne(nStateGrid))^ggamma1...
                *ttheta1*vTangebleCapitalOne(nStateGrid)^(-1/llambda1);
            RHSkt2 = vMmuTwo(nStateGrid)*aalpha2*vg2(nStateGrid)^(aalpha2*llambda2/(llambda2-1)-1)*(exp(ZI)*vLaborTwo(nStateGrid))^ggamma2...
                *ttheta2*vTangebleCapitalTwo(nStateGrid)^(-1/llambda2);
            Diff(nIndexkt2,nStateGrid,nProductivity) = RHSkt2-LHSkt2;
            
            LHSki2 = mmu1*aalpha1*vg1(nStateGrid)^(aalpha1*llambda1/(llambda1-1)-1)*(exp(ZT)*vLaborOne(nStateGrid))^ggamma1...
                *(1-ttheta1)*vIntangebleCapitalOne(nStateGrid)^(-1/llambda1);
            RHSki2 = vMmuTwo(nStateGrid)*aalpha2*vg2(nStateGrid)^(aalpha2*llambda2/(llambda2-1)-1)*(exp(ZI)*vLaborTwo(nStateGrid))^ggamma2...
                *(1-ttheta2)*vIntangebleCapitalTwo(nStateGrid)^(-1/llambda2);
            Diff(nIndexki2,nStateGrid,nProductivity) = RHSki2-LHSki2;
            
            Diff(nIndexl1,nStateGrid,nProductivity) = w*vLaborOne(nStateGrid)-mmu1*vg1(nStateGrid)^(aalpha1*llambda1/(llambda1-1))*ggamma1*(exp(ZT)*vLaborOne(nStateGrid))^ggamma1;
            
            Diff(nIndexl2,nStateGrid,nProductivity) = w*vLaborTwo(nStateGrid)-vMmuTwo(nStateGrid)*vg2(nStateGrid)^(aalpha2*llambda2/(llambda2-1))*ggamma2*(exp(ZI)*vLaborTwo(nStateGrid))^ggamma2;
            
            Diff(nIndexb,nStateGrid,nProductivity) = 1-0.04*vDebtChange(nStateGrid)+EJb/(1+ri);
            nIntangebleCapitalToday = (mSmolyakGrid(nStateGrid,2)+1)*(nIntangebleCapitalMax-nIntangebleCapitalMin)/2+nIntangebleCapitalMin;

            Diff(nIndexki,nStateGrid,nProductivity) = vg2(nStateGrid)^(aalpha2*llambda2/(llambda2-1))*(exp(ZI)*vLaborTwo(nStateGrid))^ggamma2+(1-ddeltai)*nIntangebleCapitalToday-vIntangebleCapitalTomorrow(nStateGrid);
            %Diff(nIndexmmu2,nStateGrid,nProductivity) = vMmuTwo(nStateGrid)-EJki/(1+ri);
            %{
            if (flag == 1)
                Diff(nIndexkt,nStateGrid,nProductivity) = 10^9;
            elseif (flag == 3)
                Diff(nIndexb,nStateGrid,nProductivity) = 10^9;
            end
            %}    
                
        end
        
    end
    
end

    Diff = Diff;
end

end