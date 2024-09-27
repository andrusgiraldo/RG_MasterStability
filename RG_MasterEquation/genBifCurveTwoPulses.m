%% Script to generate the two-Pulse Data
% By Andrus Giraldo and Stefan Ruschel - Last big revision 23/08/2024
%
% The FitzHugh-Nagumo profile equation is given by 
% 
%   x1'(t)= x1(t) − x1(t)^3/3 − x2(t) + I + c*(x1(t-tau1) - 2*x1(t) + x1(t-tau2))
%   x2'(t)= 0.08*(x1(t)+0.7−0.8*x2(t)) 
%
% as per http://www.scholarpedia.org/article/FitzHugh-Nagumo_model
% appended by the Master stability function
%
%   u1'(t) = A(1,1)*u1(t) + A(1,2)*u2(t) - l1*u1(t) + l2*v1(t) + c*(  cos(phi)*u1(t-tau1) + sin(phi)*v1(t-tau1) + cos(phi)*u1(t-tau2) - sin(phi)*v1(t-tau2)) 
%   v1'(t) = A(1,1)*v1(t) + A(1,2)*v2(t) - l2*u1(t) - l1*v1(t) + c*(- sin(phi)*u1(t-tau1) + cos(phi)*v1(t-tau1) + sin(phi)*u1(t-tau2) + cos(phi)*v1(t-tau2))
%   u2'(t) = A(2,1)*u1(t) + A(2,2)*u2(t) - l1*u2(t) + l2*v2(t) 
%   v2'(t) = A(2,1)*v1(t) + A(2,2)*v2(t) - l2*u2(t) - l1*v2(t) 
% where A=[1-V^2-2*c, -1; 0.08, -0.064];
%
% Parameters are (in this order) |I|,|c|,|tau1|,|tau2|,|l1|,|l2|,|phi|,|T|,|k|, ...
%
% In this script we start from a particular initial condition to find the
% the period two-pulse, and the reconstruct the bifurcation diagram as
% function of k.


%% Loading the path where the DDE-Biftool scripts are stored
clc; clear all; close all; 
addpath('../../DDE-Biftool/ddebiftool',...
    '../../DDE-Biftool/ddebiftool_extra_psol',...
    '../../DDE-Biftool/ddebiftool_extra_nmfm',...
    '../../DDE-Biftool/ddebiftool_extra_symbolic',...
    '../../DDE-Biftool/ddebiftool_utilities');
format compact

%% Set Parameter names and  function files for DDE-BifTools
parnames    =   {'I','c','tau1','tau2' 'l1', 'l2', 'phi', 'T','k'};
cind        =   [parnames;num2cell(1:length(parnames))];
ind         =   struct(cind{:});


funcs       =   set_symfuncs(@sym_fhn,'sys_tau',@()[ind.tau1, ind.tau2]);
funcs_k     =   set_symfuncs(@sym_fhn,'sys_tau',@()[ind.tau1, ind.tau2],'sys_cond', @(p)condFixK(p,ind));
funcs_ext   =   set_symfuncs(@sym_ext_fhn,'sys_tau',@()[ind.tau1, ind.tau2],'sys_cond', @(p)condTestSys(p,ind));

%% Computing the case for q=1/20

load(['./Data/SimulationMainPeriodic_inv_q_50.mat'])

plot(spatioTemporal.simul{2}.sol.y(1,:))
indices             =   [48073, 49361];
auxTwoPer.x         =   spatioTemporal.simul{2}.sol.x(indices(1):indices(2));
auxTwoPer.y(1,:)    =   spatioTemporal.simul{2}.sol.y(1,indices(1):indices(2));
auxTwoPer.y(2,:)    =   spatioTemporal.simul{2}.sol.y(2,indices(1):indices(2));
auxTwoPer.period    =   auxTwoPer.x(end)-auxTwoPer.x(1);
auxTwoPer.x         =   (auxTwoPer.x-auxTwoPer.x(1))/auxTwoPer.period;

I                   =   0.2000; 
c                   =   1.0; 
T                   =   auxTwoPer.period;
tau1                =   auxTwoPer.period/100;
tau2                =  -tau1;
k                   =   tau1/T;
l1=0; l2=0; phi=0; 

% Initial parameters;
parameters  =  cellfun(@(x)evalin('caller',x),parnames);

varargin    = {'extra_condition',true, 'matrix', 'sparse', ...
            'newton_max_iterations', 15, 'print_residual_info',1};

bounds      = {'min_bound',[],'max_step',[0,1.]};

%% Starting from a precomputed solution
sol.x       =   linspace(0,1,1000); 
sol.y(1,:)  =   interp1(auxTwoPer.x, auxTwoPer.y(1,:),sol.x,'spline');
sol.y(2,:)  =   interp1(auxTwoPer.x, auxTwoPer.y(2,:),sol.x,'spline');

dFun1       =   diff(sol.y(1,:));
dFun1       =   [dFun1,dFun1(1)];
dFun2       =   diff(sol.y(2,:));
dFun2       =   [dFun2,dFun2(1)];

nn          =   sqrt(sum(dFun1.^2+dFun2.^2)*T/1000);
dFun1       =   dFun1/nn;
dFun2       =   dFun2/nn;

sol.y(3,:)  =   dFun1;
sol.y(4,:)  =   0*sol.x;
sol.y(5,:)  =   dFun2;
sol.y(6,:)  =   0*sol.x;


sol.xe      =   [sol.x(1), sol.x(end)];
sol.ye      =   [sol.y(:,1), sol.y(:,end)];

% setting up the initial continuation run
default     =   {'degree',5,'intervals', 200,'corpar',[], 'step', 1e-2, 'indperiod',[],'stepcond',[]};
[options,pass_on]   =   dde_set_options(default,varargin,'pass_on');

if ~isfield(sol,'xe') || length(sol.xe)<2
    error('branch_from_sol:cross_section',...
        'solution structure does not contain two Poincare cross section events');
end

pmesh       =   linspace(0,1,options.degree*options.intervals+1);
profile     =   interp1(sol.x',sol.y',pmesh,'pchip')'; 
per0        =   dde_psol_create('parameter',parameters,'degree',options.degree,...
                    'profile',[],'period',T,...
                    'mesh',pmesh,'profile',profile);

per1        =   p_remesh(per0,options.degree,options.intervals);
branch      =   df_brnch(funcs_ext,[ind.l1, ind.phi, ind.l2, ind.k, ind.tau1],'psol');
branch      =   replace_branch_pars(branch,[ind.l1, ind.phi, ind.l2, ind.k, ind.tau1],pass_on);

% put in varargin by hand
branch.method.point.extra_condition         =   1;
branch.method.point.print_residual_info     =   1;
branch.method.point.newton_max_iterations   =   15;
branch.method.point.newton_nmon_iterations  =   10;
branch.method.point.matrix                  =   'sparse';

corpar      = [ind.l1, ind.l2, ind.phi, ind.tau1,ind.tau2]; 
[per2,suc]  = p_correc(funcs_ext,per1,corpar,options.stepcond,branch.method.point);
if ~suc
    error('branch_from_sol:correction',...
        '1st correction failed: suc=%d',suc);
end
branch.point    =   per2;

auxPoint                    =   branch.point;
auxPoint.parameter(ind.phi) =   branch.point.parameter(ind.phi)-options.step/1000;
[auxPoint,suc]              =   p_correc(funcs_ext,auxPoint,corpar,options.stepcond,branch.method.point);
if ~suc
    error('branch_from_sol:correction',...
        '2nd correction failed: suc=%d',suc);
end

% initialize branch
branch.point(2) = auxPoint;

[initBranch,suc]               =   SetupPsolFrom_psol(funcs_ext, branch, 2,'contpar', [ind.l1, ind.phi, ind.l2, ind.k, ind.tau1], 'extra_condition','true', 'step', 1e-3);
initBranch.parameter.max_bound =   [ind.phi, 0.01]; 
initBranch.parameter.min_bound =   [ind.phi, -0.01];
initBranch.parameter.max_step  =   [0, 1.; ind.phi, 0.01];

figure(2); clf; hold on;
initBranch = br_contn(funcs_ext,initBranch,1000);
initBranch = br_rvers(initBranch);
initBranch = br_contn(funcs_ext,initBranch,1000);
hold off;

%% Continue in wave number to obatin the main bifurcation curve shown in Figure 2 and 3 
auxbr                   =   initBranch; 
[br5,suc]               =   SetupPsolFrom_psol(funcs_ext, auxbr, 1,'contpar', [ind.tau1, ind.k, ind.tau2, ind.l1, ind.l2], 'extra_condition','true', 'step', 1e-3);
br5.parameter.min_bound =   [ind.k, 1/200];
br5.parameter.max_bound =   [ind.k, 1/12];
br5.parameter.max_step  =   [0, 1.; ind.tau1, 0.002; ind.k, 0.002];

figure(2); clf; hold on;
br5                     =   br_contn(funcs_ext,br5,2000);

BifCurveTwoPulses.branch         =   br5;
BifCurveTwoPulses.contPar        =   {'k','tau1','tau2','l1','l2'};
BifCurveTwoPulses.desp           =   'Bifurcation curve of the traveling wave with two pulses in the FH system.';

save('./Data/BifurcationCurveTwoPulses.mat', 'BifCurveTwoPulses');