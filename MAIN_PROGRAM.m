% --- SOMA Simple Program --- Version: SOMA T3A (V1.0) August 25, 2020 ----
% ------ Written by: Quoc Bao DIEP ---  Email: diepquocbao@gmail.com   ----
% -----------  See more details at the end of this file  ------------------
clearvars; %clc;
disp('Hello! SOMA T3A is working, please wait... ')
tic                                                                         % Start the timer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  U S E R    D E F I N I T I O N
% Define the Cost function, if User wants to optimize another function,
% Please change the function name, for example:  CostFunction = @(pop)      schwefelfcn(pop);
%                                           or:  CostFunction = @(pop)      ackleyfcn(pop);
%                                           or:  CostFunction = @(pop)      periodicfcn(pop);
%CostFunction = @(pop)     schwefelfcn(pop);                                % use this line if: pop = PopSize x Dimension
CostFunction  = @(pop)     schwefelfcn(pop');                               % use this line if: pop = Dimension x PopSize
    %--------------- Initial Control Parameters of SOMA -------------------
            SOMApara.PopSize    = 100;                                      % The population size of SOMA
            SOMApara.N_jump     = 45;                                       % The number of jump of each individual
			SOMApara.m     		= 10;                                       % The parameter m
			SOMApara.n     		= 5;                                        % The parameter n
			SOMApara.k     		= 10;                                       % The parameter k
    %----------------------------------------------------------------------
            Info.dimension      = 10;                                       % Define the dimension of the problem
            Info.Search_Range   = [-500 500];                               % Define the search space (lower-upper)
            Info.Max_FEs        = 1e4*Info.dimension;      	                % The unique stop condition of the algorithm
    %       E N D    O F     U S E R
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [Best , FEs , Mig]      = SOMA_T3A_v1(Info,SOMApara,CostFunction);  % Call SOMA T3A Algorithm
time = toc;                                                                 % Stop the timer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Show the information to User
disp(['Stop at Migration :  ' num2str(Mig)])
disp(['The number of FEs :  ' num2str(FEs)])
disp(['Processing time   :  ' num2str(time)])
disp(['The best cost     :  ' num2str(Best.Value)])
disp(['Solution values   :  ']), disp(Best.Positon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This algorithm is programmed according to the descriptions in the papers 
% listed below:

% Link of paper: https://ieeexplore.ieee.org/abstract/document/8790202/
% Diep, Q.B., 2019, June. Self-Organizing Migrating Algorithm Team To Team 
% Adaptive???SOMA T3A. In 2019 IEEE Congress on Evolutionary Computation (CEC)
% (pp. 1182-1187). IEEE.

% The control parameters PopSize, N_jump, m, n, and k are closely related 
% and greatly affect the performance of the algorithm. Please refer to the 
% above paper to use the correct control parameters.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LIST OF COST FUNCTIONS (Need to optimize):

% You will replace these functions with your own that describes the problem you need to optimize for.
%--------------------------------------------------------------------------
% Schwefel function:
% schwefelfcn accepts a matrix of size M-by-N and returns a vetor SCORES
% of size M-by-1 in which each row contains the function value for each row of X.
function scores = schwefelfcn(x)
    n = size(x, 2);
    scores = 418.9829 * n - (sum(x .* sin(sqrt(abs(x))), 2));
end
%--------------------------------------------------------------------------
% Ackley function:
% ackleyfcn accepts a matrix of size M-by-N and returns a vetor SCORES
% of size M-by-1 in which each row contains the function value for each row of X.
function scores = ackleyfcn(x)
    n = size(x, 2);
    ninverse = 1 / n;
    sum1 = sum(x .^ 2, 2);
    sum2 = sum(cos(2 * pi * x), 2);
    scores = 20+exp(1)-(20*exp(-0.2*sqrt(ninverse*sum1)))-exp(ninverse*sum2);
end
%--------------------------------------------------------------------------
% Periodic function:
% periodicfcn accepts a matrix of size M-by-N and returns a vetor SCORES
% of size M-by-1 in which each row contains the function value for each row of X.
function scores = periodicfcn(x)
    sin2x = sin(x) .^ 2;
    sumx2 = sum(x .^2, 2);
    scores = 1 + sum(sin2x, 2) -0.1 * exp(-sumx2);
end
