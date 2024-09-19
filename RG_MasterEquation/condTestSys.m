%% Definition of the zero condition in DDE-BifTools format
% By Andrus Giraldo and Stefan Ruschel - Last big revision 23/08/2024

function [res,J]=condTestSys(p, ind)
%CONDTESTSYS Boundary and zero conditions for the Master Stability Problem
%
%   Input Variables:
%       p:      DDEBifTool branch point solution structure.
%       ind:    Structure with indexes of the parameters. 

% Number of extra conditions is m
m = 4; res = zeros(m,1); J = repmat(p_axpy(0,p,[]),m,1);

% Get information at t=0
[x,Jx]          =   dde_coll_eva(p.profile,p.mesh,0,p.degree,'kron',true);

% Fix L2-norm of the components 3-6 to be one
[auxM, ~]       =   size(p.profile); 
auxWeight       =   repmat([diff(p.mesh),0], [auxM,1]) * p.period;
auxEval         =   p.profile.*auxWeight;
C1              =   diag(repmat([0 0 1 1 1 1],1,length(p.mesh)));
res(1)          =   auxEval(:)'*C1*p.profile(:)-1;
J(1).profile    =   reshape(2*auxEval(:)'*C1,size(p.profile));  
J(1).period     =   (res(1)+1)/p.period;

% Fix value of the fourth component to zero at zero
C2              =   [0 0 0 1 0 0];
res(2)          =   C2*x;
J(2).profile    =   reshape(C2*Jx,size(p.profile));

% Fix the values of tau1 and tau2 to be equal in absolute value
res(3)                      =   p.parameter(ind.tau1) + p.parameter(ind.tau2);
J(3).parameter(ind.tau1)    =   1;
J(3).parameter(ind.tau2)    =   1;

% Define the relation tau/period = k
res(4)                      =   p.parameter(ind.tau1) - p.parameter(ind.k)*p.period;
J(4).parameter(ind.tau1)    =   1;
J(4).parameter(ind.k)       =   -p.period;
J(4).period                 =   -p.parameter(ind.k);