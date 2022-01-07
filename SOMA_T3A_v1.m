% --- SOMA Simple Program --- Version: SOMA T3A (V1.0) August 25, 2020 ----
% ------ Written by: Quoc Bao DIEP ---  Email: diepquocbao@gmail.com   ----
% -----------  See more details at the end of this file  ------------------
function [Best , FEs , Mig] = SOMA_T3A_v1(Info , SOMApara , CostFunction)
    % -------------- Extract Information ----------------------------------
    Max_FEs           = Info.Max_FEs;
    dimension         = Info.dimension;
    % -------------- The search range of the function ---------------------
    VarMin            = Info.Search_Range(1);
    VarMax            = Info.Search_Range(2);
    % -------------- Initial Parameters of SOMA ---------------------------
    PopSize           = SOMApara.PopSize;
    N_jump            = SOMApara.N_jump;
    m                 = SOMApara.m;
    n                 = SOMApara.n;
    k                 = SOMApara.k;
    % --------------------- Create Initial Population ---------------------
    pop               = VarMin + rand(dimension,PopSize)*(VarMax - VarMin);
    fit               = CostFunction(pop);
    FEs               = PopSize;
    [global_cost, id] = min(fit);
    global_pos        = pop(:,id);
    % ---------------- SOMA MIGRATIONS ------------------------------------
    Mig         = 0;
    while (FEs+N_jump < Max_FEs)
        Mig     = Mig  + 1;
		% ------------ Update PRT and Step parameters ---------------------
        PRT     = 0.05 + 0.90*(FEs/Max_FEs);
        Step    = 0.15 - 0.08*(FEs/Max_FEs);
        % ------------ Migrant selection: m -------------------------------
        M       = randperm(PopSize,m);
        sub_M   = [fit(M)' ; pop(:,M) ; M];
        sub_M   = sortrows(sub_M')';
        fit_M   = sub_M(1,:);
        pop_M   = sub_M(2:end-1,:);
        order_M = sub_M(end,:);
        % ------------ movement of n Migrants -----------------------------
        for j  = 1 : n
            Migrant = pop_M(:,j);
            % ------------ Leader selection: k ----------------------------
            K       = randperm(PopSize,k);
            sub_K   = [fit(K)' ; pop(:,K) ; K];
            sub_K   = sortrows(sub_K')';
            pop_K   = sub_K(2:end-1,:);
            order_K = sub_K(end,:);
            Leader  = pop_K(:,1);
            % ------------- Moving process --------------------------------
            coincide = order_M(j) ~= order_K(1);
            if  coincide
                offs_path     = [];
                for move      = 1 : N_jump
                    nstep     = move * Step;
                    %----- SOMA Mutation ----------------------------------
                    PRTVector = rand(dimension,1) < PRT;
                    %----- SOMA Crossover ---------------------------------
                    offspring =  Migrant + (Leader - Migrant)*nstep.*PRTVector;
                    offs_path = [offs_path   offspring];
                end % END JUMPING
                %----- Checking Boundary and Replaced Outsize Individuals -
                number_offs_path = size(offs_path,2);
                for cl = 1 : number_offs_path
                    for rw = 1 : dimension
                        if  (offs_path(rw,cl) < VarMin) || (offs_path(rw,cl) > VarMax)
                             offs_path(rw,cl) = VarMin  +  rand*(VarMax - VarMin);
                        end
                    end
                end
                %----- SOMA Re-Evaluate Fitness Fuction -------------------
                new_cost               = CostFunction(offs_path);
                FEs                    = FEs + number_offs_path;
                %----- SOMA Accepting: Place Best Individual to Population-
                [min_new_cost, idz]    = min(new_cost);
                if  min_new_cost      <= fit_M(j)
                    indi_replace       = offs_path(:,idz);
                    pop(:,order_M(j))  = indi_replace;
                    fit(order_M(j))    = min_new_cost;
                    %----- SOMA Update Global_Leader ----------------------
                    if  min_new_cost   < global_cost
                        global_cost    = min_new_cost;
                        global_pos     = indi_replace;
                    end
                end
            end % END indi_moving ~= Target
        end  % END PopSize   (For Loop)
    end   % END MIGRATIONS (While Loop)
    Best.Value      = global_cost;
    Best.Positon    = global_pos;
end
%--------------------------------------------------------------------------
% This algorithm is programmed according to the descriptions in the papers 
% listed below:

% Link of paper: https://ieeexplore.ieee.org/abstract/document/8790202/
% Diep, Q.B., 2019, June. Self-Organizing Migrating Algorithm Team To Team 
% Adaptive–SOMA T3A. In 2019 IEEE Congress on Evolutionary Computation (CEC)
% (pp. 1182-1187). IEEE.

% The control parameters PopSize, N_jump, m, n, and k are closely related 
% and greatly affect the performance of the algorithm. Please refer to the 
% above paper to use the correct control parameters.