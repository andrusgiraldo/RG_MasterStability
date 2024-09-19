%% Simulation of the bi-direction FH system
% By Andrus Giraldo and Stefan Ruschel - Last big revision 23/08/2024

function [solAux] = compSimulationBiDi(auxPerPoint,M,tf,shift,scaledPerturb,speCondition)
%COMPSIMULATIONBIDI Simulates network of FH nodes that are coupled
%   throught their nearest neighbour. This is in a ring network 
%
%   As a starting initial condition, we take a periodic solution that was
%   used for our Master Stability curve computations, and use it to start 
%   each node accordingly. Additionaly, from this periodic solution, we 
%   construct the corresponding network size that supports it.
%
%   Input Variables:
%       auxPerPoint:    Function with the linearized internal dynamics 
%                       of each node. This variable has as input "t"
%       M:              Number of copies
%       tf:             Final time of integration
%       shift:          Quantity to shit the starting time along the
%                       periodic solution, where the initial conditions 
%                       will be sampled
%       scaledPerturb:  Perturbation in the magnitude of the initial cond.
%       speCondition:   If speCondition is non-empty, then the initial
%                       condition for simulations are going to be the 
%                       entries of this vector.

period      =   auxPerPoint.period;
I           =   auxPerPoint.parameter(1);
c           =   auxPerPoint.parameter(2);
tau         =   auxPerPoint.parameter(3);

auxData     =   [auxPerPoint.mesh', auxPerPoint.profile(1,:)', auxPerPoint.profile(2,:)'];

% Defining the single unit fields and the coupling of the filds
fSingleField   = @(U) [U(1)-U(1)^3/3-U(2)+I; 0.08*(U(1)+0.7-0.8*U(2))];
gCoupling      = @(U) [c*U(1);0];


% Creating the one-directional coupling matrix 
nDimUnit            =   2;
nDimNetwork         =   M*round(period/tau);
matrixCoupling      =   zeros(nDimNetwork,nDimNetwork);
for i=1:nDimNetwork
    matrixCoupling(mod(i-1,nDimNetwork)+1,mod(i,nDimNetwork)+1)     =   1;
    matrixCoupling(mod(i-1,nDimNetwork)+1,mod(i-2,nDimNetwork)+1)   =   1;
    matrixCoupling(i,i)                                             =   -2;
end

fFieldNet   = @(t,U) createNetworkDynamics(fSingleField,gCoupling,matrixCoupling,nDimUnit,U);

% Setup the initial conditions
if isempty(speCondition)==1
    % Creating a Matlab interpolating object used for the initial conditions
    perFun1     = pchip(auxData(:,1), auxData(:,2));
    perFun2     = pchip(auxData(:,1), auxData(:,3));
    initCond    = [zeros(1,nDimNetwork);zeros(1,nDimNetwork)];
    % Sample the initial conditions using the periodic solution
    for i=1:nDimNetwork
        if isempty(scaledPerturb)
            initCond(1,i) = ppval(perFun1,mod(-(i-1)/nDimNetwork*M  + shift,1));
            initCond(2,i) = ppval(perFun2,mod(-(i-1)/nDimNetwork*M  + shift,1));
        else
            initCond(1,i) = ppval(perFun1,mod(-(i-1)/nDimNetwork*M  + shift,1))*(scaledPerturb(i));
            initCond(2,i) = ppval(perFun2,mod(-(i-1)/nDimNetwork*M  + shift,1));
        end
    end
else
    initCond    = speCondition;
end
options1    =   odeset('RelTol',1e-08,'AbsTol',1e-07*ones(1,nDimNetwork*nDimUnit)); 
solAux      =   ode45(fFieldNet,[0,tf],initCond(:),options1);  
end

