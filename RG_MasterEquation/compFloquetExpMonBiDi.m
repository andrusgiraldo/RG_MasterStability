%% Computation of Floquet Exponent by computing the Monondromy Matrix
% By Andrus Giraldo and Stefan Ruschel - Last big revision 23/08/2024

function [floqM, floqE, auxVec] = compFloquetExpMonBiDi(auxPerPoint,M)
%COMPFLOQUETEXPMOBIDI Screates a fFlowNetwork variable where every nDimUnit
%   rows has the linearized dynamics of a node.
%
%   We are considering the following Vector field
%       w_i'= A_i(t) *w_i + sum_j A_{ij}*B_j(t)*w_j, w_i \in R^nDimUnit
%   This function output corresponds to the value of the right-hand side of
%   the equation above. Here, A_i and B_j are T-periodic matrices (the equation 
%   above is a Floquet problem).
%
%   Furthermore, we are considering a periodic traveling wave solution;
%   hence, the evaluation in time of matrices of A_i and B_j will depend on
%   two matrices A and B that are evaluated at different fractional 
%   proportions of the period from a reference time t0.
%
%   For example, one could have the following for w_1 at time t=t0
%       w_1'= A(t0+k_1*T) *w_1 + sum_j A_{ij}*B(t0+T*k_j)*w_j
%   where, k_j is a number between 0 and 1. Notice that there is only one 
%   matrix A and one matrix B instead of matrices A_i and B_j. The matrices 
%   are evaluated at time k_j*T.
%
%   Input Variables:
%       auxPerPoint:    Function with the linearized internal dynamics 
%                       of each node. This variable has as input "t"
%       M:              Number of copies


period      =   auxPerPoint.period;
I           =   auxPerPoint.parameter(1);
c           =   auxPerPoint.parameter(2);
tau         =   auxPerPoint.parameter(3);

auxData     =   [auxPerPoint.mesh', auxPerPoint.profile(1,:)', auxPerPoint.profile(2,:)'];

% Creating Spline condition
perFun      =   pchip(auxData(:,1)*period, auxData(:,2)); % We only need the first component
perFunEval  =   @(t) ppval(perFun,t); 


% Definition of the linear couplings
fLinearSingleField  = @(t) [1-perFunEval(t)^2, - 1; 0.08, -0.064];                         
gLinearCoupling     = @(t) [c, 0; 0, 0];

% Creating the one-directional coupling matrix 
nDimUnit            =   2;
nDimNetwork         =   M*round(period/tau);
matrixCoupling      =   zeros(nDimNetwork,nDimNetwork);
for i=1:nDimNetwork
    matrixCoupling(mod(i-1,nDimNetwork)+1,mod(i,nDimNetwork)+1)     =   1;
    matrixCoupling(mod(i-1,nDimNetwork)+1,mod(i-2,nDimNetwork)+1)   =   1;
    matrixCoupling(i,i)                                             =   -2;
end

auxIndTimes         =   linspace(0,1,round(period/tau)+1);
indTimes            =   repmat(-1*auxIndTimes(1:end-1),1,M);


fLinearFieldNet     =   @(t,U) createLinearizedNetworkDynamics(fLinearSingleField,gLinearCoupling,matrixCoupling,nDimUnit,U,indTimes,t,period);
options1            =   odeset('RelTol',1e-8,'AbsTol',1e-7*ones(1,nDimNetwork*nDimUnit)); 

initCond            =   eye(nDimUnit*nDimNetwork);
monMatrix           =   zeros(nDimUnit*nDimNetwork);
for i = 1:nDimUnit*nDimNetwork
    disp(['Basis vector ',num2str(i),' of ',num2str(nDimUnit*nDimNetwork)]);
    solAux              =   ode45(fLinearFieldNet,[0,period],initCond(:,i),options1); %ode45
    monMatrix(:,i)      =   solAux.y(:,end);  
end

[auxVec, auxEig]    =   eig(monMatrix,"balance");
floqM               =   diag(auxEig);
floqE               =   log(floqM)/period;
end

