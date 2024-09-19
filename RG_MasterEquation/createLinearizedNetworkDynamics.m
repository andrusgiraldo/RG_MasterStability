%% Definition of the linear stability problem of a periodic solution
% By Andrus Giraldo and Stefan Ruschel - Last big revision 23/08/2024

function [fFlowNetwork] = createLinearizedNetworkDynamics(fLinearFunctionUnit,gLinearCoupling,matrixCoupling,nDimUnit,U,indTimes,time,period)
%CREATENETWORKDYNAMICS creates a fFlowNetwork variable where every nDimUnit
%   rows has the linearized dynamics of node.
%
%   We are considering the following Vector field
%       w_i'= A_i(t) *w_i + sum_j A_{ij}*B_j(t)*w_j, w_i \in R^nDimUnit
%   This function output corresponds to the value of the right-hand side of
%   the equation above. Here, A_i and B_j are T-periodic matrices (the equation 
%   above is a Floquet  problem).
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
%       fLinearFunctionUnit:  Function with the linearized internal dynamics 
%                             of each node. This variable has as input "t"
%                             and output the matrix --A(t)--
%       gCoupling:            Function with the linearized coupling dynamics 
%                             of each node. This variable has as input "t"
%                             and output the matrix --B(t)--
%       matrixCoupling:       Adjancency matrix of the network --A_{ij}--.
%       nDimUnit:             Dimension of the internal dynamics --dim(u_j)--.
%       U:                    Vector with the values --w--.
%       indTimes:             Vector with values between 0 and 1 --k--.
%       time:                 Time where the dynamics are evaluated (t=t0).
%       period:               The minimal period of matrices A and B.

[nDimNetwork,~]   = size(matrixCoupling);
fFlowNetwork      = zeros(nDimUnit*nDimNetwork,1);
for i=1:nDimNetwork
    fFlowNetwork((i-1)*nDimUnit+1:i*nDimUnit) = fLinearFunctionUnit(mod(time+period*indTimes(i),period))*U((i-1)*nDimUnit+1:i*nDimUnit);
    for j=1:nDimNetwork
        if (matrixCoupling(i,j)~=0)
            fFlowNetwork((i-1)*nDimUnit+1:i*nDimUnit) = fFlowNetwork((i-1)*nDimUnit+1:i*nDimUnit) + matrixCoupling(i,j)*gLinearCoupling(mod(time+period*indTimes(j),period))*U((j-1)*nDimUnit+1:j*nDimUnit);
        end
    end
end