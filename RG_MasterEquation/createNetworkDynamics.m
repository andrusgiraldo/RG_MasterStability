%% Definition of the System of ODEs in a given network
% By Andrus Giraldo and Stefan Ruschel - Last big revision 23/08/2024

function [fFlowNetwork] = createNetworkDynamics(fFunctionUnit,gCoupling,matrixCoupling,nDimUnit,U)
%CREATENETWORKDYNAMICS creates a fFlowNetwork variable where every nDimUnit
%   rows have the value of the derivative in time of a node.
%
%   We are considering that the network is written in vector field form
%       u_i'=f(u_i)+ sum_j A_{ij}*H(u_i,u_j), u_i \in R^nDimUnit
%   This function output corresponds to the value of the right-hand side of
%   the equation above.
%
%   Input Variables:
%       fFunctionUnit:  Function with the internal dynamics of each node --f--.
%       gCoupling:      Function with the functional relation of a node and 
%                       its neighbourgh --H(x_i,x_j)--.
%       matrixCoupling: Adjancency matrix of the network --A_{ij}--.
%       nDimUnit:       Dimension of the internal dynamics --dim(u_j)--.
%       U:              Vector with the value of the nodes of the network 
%                       at a time t  --(u_1,u_2,...,u_n)--.

[nDimNetwork,~]   = size(matrixCoupling);
fFlowNetwork      = zeros(nDimUnit*nDimNetwork,1);

for i=1:nDimNetwork
    fFlowNetwork((i-1)*nDimUnit+1:i*nDimUnit) = fFunctionUnit(U((i-1)*nDimUnit+1:i*nDimUnit));
    for j=1:nDimNetwork
        if (matrixCoupling(i,j)~=0)
            fFlowNetwork((i-1)*nDimUnit+1:i*nDimUnit) = fFlowNetwork((i-1)*nDimUnit+1:i*nDimUnit) + matrixCoupling(i,j)*gCoupling(U((j-1)*nDimUnit+1:j*nDimUnit));
        end
    end
end