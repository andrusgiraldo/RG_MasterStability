%% Script to generate the main bifurcation curve shown in Figure 2 and 3 
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
% We first find a periodic solution when varying tau1 while
% tau2=0. Then, we perform continuation in tau2 to make both delay and advanced
% variable equal in modulus. We then perform numerical continuation to initialize our
% Master Equation Curve and continue in the wave number k. See below for more details.

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


%% Initialization of branch of non-trivial equilibria
% This branch of non-trivial equilibria will exhibit a Hopf bifurcation,
% which we will use to compute an initial periodic solution.

% Initialization of parameter and initial conditions. Notice that
% tau1=tau2=1

I=0.8750; c=1; tau1=1; tau2=1; l1=0; l2=0; phi=0; T=1; k=1;
x0          =   [0;0.8750];
par0        =   cellfun(@(x)evalin('caller',x),parnames);

contpar     =   ind.I;
bounds      =   {'min_bound',[],'max_step',[0,1.]};
nontriv_eqs =   SetupStst(funcs,'x',x0,'parameter',par0,'step', 0.1,...
                'contpar',contpar,'max_step',[contpar,0.3],bounds{:});

% Compute, find stability and bifurcations of non-trivial equilibria as the
% parameter -I- is varied
disp('Trivial equilibria');
figure(1);clf;ax1=gca;
nontriv_eqs                         =   br_contn(funcs,nontriv_eqs,20,'plotaxis',ax1);
[nontriv_eqs,~,ind_hopf,bif1types]  =   LocateSpecialPoints(funcs,nontriv_eqs);
nunst_eqs                           =   GetStability(nontriv_eqs);

% Find index -ind_hopf- in the solution branch where a Hopf bifurcation occurs
ind_hopf                            =   find(diff(nunst_eqs)~=0);
fprintf('Hopf bifurcation near point %d\n', ind_hopf);

%% Branch off at degenerate Hopf bifurcation
% That is, we compute a periodic solution branch as parameter I varies.

disp('Branch off at Hopf bifurcation');
fprintf('Initial correction of periodic orbits at Hopf:\n');

[per_orb,suc]   =   SetupPsol(funcs,nontriv_eqs,ind_hopf(1),'step', 1e-3,...
                        'print_residual_info',1,'intervals',100,'degree',10,...
                        'min_bound',[contpar,0.4],'max_step',[0, 1.1],...
                        'matrix','sparse', 'eigmatrix', 'sparse');
figure(1); clf; hold on;
per_orb         =   br_contn(funcs,per_orb,1000);
hold off;

%% We then continue in tau2 until tau2 = -1.
[br,suc]                =   SetupPsolFrom_psol(funcs,per_orb,length(per_orb.point),'contpar', ind.tau2, 'step', 1e-3);
br.parameter.min_bound  =   [ind.tau2, -1];
br.parameter.max_step   =   [0, 1.];

figure(2); clf; hold on;
br                      = br_rvers(br);
br                      = br_contn(funcs,br,100);
hold off;

%% We then continue the periodic solution under the additional conditions 
% tau1=-tau2, and the relation tau/period = k

% Here, we initialize parameter k
auxbr                               =   br;
auxbr.point(end).parameter(ind.k)   =   br.point(end).parameter(tau1)/br.point(end).period;

[br1,suc]                           =   SetupPsolFrom_psol(funcs_k, br, length(br.point),'contpar', [ind.I, ind.k, ind.tau1], 'extra_condition','true');
br1.parameter.min_bound             =   [ind.I, 0.2];
br1.parameter.max_step              =   [0, 1.];

figure(2); clf; hold on;
br1                                 =   br_rvers(br1);
br1                                 =   br_contn(funcs_k,br1,200);
hold off;

%% Do a continuation run along k while tau1 and tau2 are vary
[br2,suc]                   =   SetupPsolFrom_psol(funcs_k, br1, length(br1.point),'contpar', [ind.tau1, ind.k, ind.tau2], 'extra_condition','true');
br2.parameter.min_bound     =   [ind.k, 1/50];
br2.parameter.max_bound     =   [ind.k, 1/20];
br2.parameter.max_step      =   [0, 0.1];

figure(2); clf; hold on;
br2                         =   br_contn(funcs_k,br2,1000);
delnum                      =   length(br2.point)-1;

br2.parameter.min_bound     =   [ind.k, 1/20];
br2.parameter.max_bound     =   [];
br2                         =   br_contn(funcs_k,br2,1000);
br2.point(1:delnum)         =   [];
hold off;


%% Set up the extended system for the computation of the master stability 
% equation and compute trivial stability curve

loc = 2;
auxbr                       =   br2;

% Constructing the trivial bundle
pt                          =   br2.point(loc);
dFun1                       =   diff(interp1(pt.mesh, pt.profile(1,:),linspace(0,1,1000),'pchip'));
dFun1                       =   [dFun1,dFun1(1)];
dFun2                       =   diff(interp1(pt.mesh, pt.profile(2,:),linspace(0,1,1000),'pchip'));
dFun2                       =   [dFun2,dFun2(1)];
nn                          =   sqrt(sum(dFun1.^2+dFun2.^2)*pt.period/1000);
dFun1                       =   dFun1/nn;
dFun2                       =   dFun2/nn;

% Appending the trivial bundle to the initial point for continuation
pt.profile(3,:)             =   interp1(linspace(0,1,1000), dFun1, pt.mesh,'pchip');
pt.profile(4,:)             =   0*pt.mesh;
pt.profile(5,:)             =   interp1(linspace(0,1,1000), dFun2, pt.mesh,'pchip');
pt.profile(6,:)             =   0*pt.mesh;
auxbr.point(loc)            =   pt; 

% Continue the extended system to phi=|0.01|
[br3,suc]                   =   SetupPsolFrom_psol(funcs_ext, auxbr, loc,'contpar', [ind.l1, ind.phi, ind.l2, ind.k, ind.tau1], 'extra_condition','true', 'step', 1e-4);
br3.parameter.max_bound     =   [ind.phi, 0.01]; 
br3.parameter.min_bound     =   [ind.phi, -0.01];
br3.parameter.max_step      =   [0, 0.1];

figure(2); clf; hold on;
br3                         =   br_rvers(br3);
br3                         =   br_contn(funcs_ext,br3,1000);
br3                         =   br_rvers(br3);
br3                         =   br_contn(funcs_ext,br3,1000);
hold off;

%% Continue in wave number to obatin the main bifurcation curve shown in Figure 2 and 3 
auxbr                   =   br3; 
[br5,suc]               =   SetupPsolFrom_psol(funcs_ext, auxbr, 1,'contpar', [ind.tau1, ind.k, ind.tau2, ind.l1, ind.l2], 'extra_condition','true', 'step', 1e-3);
br5.parameter.min_bound =   [ind.k, 1/200];
br5.parameter.max_bound =   [ind.k, 1/12];
br5.parameter.max_step  =   [0, 1.; ind.tau1, 0.002; ind.k, 0.002];

figure(2); clf; hold on;
br5                     =   br_contn(funcs_ext,br5,2000);
br5                     =   br_rvers(br5);
br5                     =   br_contn(funcs_ext,br5,2000);

BifCurve.branch         =   br5;
BifCurve.contPar        =   {'k','tau1','tau2','l1','l2'};
BifCurve.desp           =   'Bifurcation curve of the traveling wave in the FH system.';

save('./Data/BifurcationCurve.mat', 'BifCurve');