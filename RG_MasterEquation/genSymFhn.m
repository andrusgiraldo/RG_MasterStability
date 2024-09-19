%% Create right-hand side and derivatives from symbolic expressions
% By Andrus Giraldo and Stefan Ruschel - Last big revision 23/08/2024

% This script generates two function files -- sym_fhn & sym_ext_fhn
% corresponding to the FitzHugh-Nagumo profile equation (with addtional
% advanced and delayed feedback terms) (sym_fhn), and the extension
% with the Master Stability Equation appended (sym_ext_fhn). 

% The FitzHugh-Nagumo profile equation is given by 
% 
% $$ x1'(t)= x1(t) − x1(t)^3/3 − x2(t) + I + c*(x1(t-tau1) - 2*x1(t) + x1(t-tau2))$$ 
% $$ x2'= 0.08*(x1(t)+0.7−0.8*x2(t)) $$
%
% as per http://www.scholarpedia.org/article/FitzHugh-Nagumo_model
%
% The profile equation can be written in the form X=(x1,y2)
% X'(t) = F(X(t)) + C*X(t-tau1) + C*X(t-tau1)
% where DF(X(t))=[1-V^2-2*c, -1; 0.08, -0.064]; C=[c,0;0,0]

% We append Master Stability Equation (Linearisation of PE with additional parameters l1, l2, phi, sigma)
% Y'(t) = A(t)Y(t) - (l1+il2)*Y(t) + exp(-i*phi)*C*Y(t - tau1) + exp(i*phi)*C*Y(t - tau1)
% where A(t)=DF(X(t)) is Jacobian of PE wrt to instantaneous argument.
%
% Y=U+i*V satisfy 
% U'(t) = A(t)*U(t) - l1*U(t) + l2*V(t) + C*(  cos(phi)*U(t-tau1) + sin(phi)*V(t-tau1)) + C*(cos(phi)*U(t-tau2) - sin(phi)*V(t-tau2))
% V'(t) = A(t)*V(t) - l2*U(t) - l1*V(t) + C*(- sin(phi)*U(t-tau1) + cos(phi)*V(t-tau1)) + C*(sin(phi)*U(t-tau2) + cos(phi)*V(t-tau2))
% where U=(u1,u2), V=(v1,v2) are the real parts, and imaginary parts of the
% first and second component of Y respecitvely

% We have 
% u1'(t) = A1*u1(t) + A2*u2(t) - l1*T*u1(t) + l2*T*v1(t) + c*(  cos(phi)*u1(t-tau1) + sin(phi)*v1(t-tau1) + cos(phi)*u1(t-tau2) - sin(phi)*v1(t-tau2)) 
% v1'(t) = A1*v1(t) + A2*v2(t) - l2*T*u1(t) - l1*T*v1(t) + c*(- sin(phi)*u1(t-tau1) + cos(phi)*v1(t-tau1) + sin(phi)*u1(t-tau2) + cos(phi)*v1(t-tau2))
% u2'(t) = A3*u1(t) + A4*u2(t) - l1*T*u2(t) + l2*T*v2(t) 
% v2'(t) = A3*v1(t) + A4*v2(t) - l2*T*u2(t) - l1*T*v2(t) 
%
% Parameters are (in this order) |I|,|c|,|tau1|,|tau2|,|l1|,|l2|,|phi|,|T|,|k|, ...

%% Loading the path where the DDE-Biftool scripts are stored
addpath('../../../DDE-Biftool/ddebiftool',...
    '../../../DDE-Biftool/ddebiftool_extra_psol',...
    '../../../DDE-Biftool/ddebiftool_extra_nmfm',...
    '../../../DDE-Biftool/ddebiftool_extra_symbolic',...
    '../../../DDE-Biftool/ddebiftool_utilities');

%% Set number of delays and parameter names
ntau     =  2;
parnames =  {'I','c','tau1','tau2' 'l1', 'l2', 'phi', 'T','k'};
cind     =  [parnames;num2cell(1:length(parnames))];
ind      =  struct(cind{:});

%% Define system using symbolic algebra using arbitrary variable names
%  Here, we are defining first the advanced-delayed differential equation
%  together with its Master Stability equation.

syms(parnames{:});
x        =  sym('x',[6,ntau+1]);
par      =  sym(parnames);

A        = [1-x(1,1)^2-2*c, -1; 0.08, -0.064];
f        = [x(1,1) - x(1,1)^3/3 - x(2,1) + I + c*(x(1,2) - 2*x(1,1) + x(1,3));
            0.08*(x(1,1) + 0.7 - 0.8*x(2,1));
            A(1,1)*x(3,1) + A(1,2)*x(5,1) - l1*x(3,1) + l2*x(4,1) + c*(  cos(phi)*x(3,2) + sin(phi)*x(4,2) + cos(phi)*x(3,3) - sin(phi)*x(4,3));  
            A(1,1)*x(4,1) + A(1,2)*x(6,1) - l2*x(3,1) - l1*x(4,1) + c*(- sin(phi)*x(3,2) + cos(phi)*x(4,2) + sin(phi)*x(3,3) + cos(phi)*x(4,3));
            A(2,1)*x(3,1) + A(2,2)*x(5,1) - l1*x(5,1) + l2*x(6,1); 
            A(2,1)*x(4,1) + A(2,2)*x(6,1) - l2*x(5,1) - l1*x(6,1)]; 
     
%% Differentiate and generate code
% Here, we export just the advanced-delayed FitzHugh-Nagumo equation to sym_fhn.m
[~, ~]   = dde_sym2funcs(f(1:2),x(1:2,:),par,'filename','sym_fhn','directional_derivative',false);

% Here, we export just the advanced-delayed FitzHugh-Nagumo equation to sym_ext_fhn.m
[~, ~]   = dde_sym2funcs(f,x,par,'filename','sym_ext_fhn','directional_derivative',false);