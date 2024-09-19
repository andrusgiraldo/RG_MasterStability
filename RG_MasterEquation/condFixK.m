%% Definition of tau1=-tau2 and the k zero condition in DDE-BifTools format
% By Andrus Giraldo and Stefan Ruschel - Last big revision 23/08/2024

function [res,J] = condFixK(p, ind)
%CONDFIXQ Zero condition that ensure that both tau1 and tau2 are equal to
%   zero in absolute value, and that k = tau1/period
%
%   Input Variables:
%       p:      DDEBifTool branch point solution structure.
%       ind:    Structure with indexes of the parameters. 

% Number of extra conditions is m
m=2; res = zeros(m,1); J = repmat(p_axpy(0,p,[]),m,1);

% Fix the values of tau1 and tau2 to be equal in absolute value
res(1)                      =   p.parameter(ind.tau1) + p.parameter(ind.tau2);
J(1).parameter(ind.tau1)    =   1;
J(1).parameter(ind.tau2)    =   1;

% Define the relation tau/period = k
res(2)                      =   p.parameter(ind.tau1) - p.parameter(ind.k)*p.period;
J(2).parameter(ind.tau1)    =   1;
J(2).parameter(ind.k)       =   -p.period;
J(2).period                 =   -p.parameter(ind.k);