%% Script to generate the main the Master stability curves shown in Figures 2 and 3 
% By Andrus Giraldo and Stefan Ruschel - Last big revision 23/08/2024
%
% The FitzHugh-Nagumo ring equation is given by 
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
% In this script we load the bifurcation curve file, and then select the
% curve from which the master stability curve are going to be computed. 
% See below for more details.

%% Loading the path where the DDE-Biftool scripts are stored
clc; clear all; close all; 
addpath('../../../DDE-Biftool/ddebiftool',...
    '../../../DDE-Biftool/ddebiftool_extra_psol',...
    '../../../DDE-Biftool/ddebiftool_extra_nmfm',...
    '../../../DDE-Biftool/ddebiftool_extra_symbolic',...
    '../../../DDE-Biftool/ddebiftool_utilities');
format compact

%% Set Parameter names and  function files for DDE-BifTools
parnames    =   {'I','c','tau1','tau2' 'l1', 'l2', 'phi', 'T','k'};
cind        =   [parnames;num2cell(1:length(parnames))];
ind         =   struct(cind{:});

funcs       =   set_symfuncs(@sym_fhn,'sys_tau',@()[ind.tau1, ind.tau2]);
funcs_k     =   set_symfuncs(@sym_fhn,'sys_tau',@()[ind.tau1, ind.tau2],'sys_cond', @(p)condFixK(p,ind));
funcs_ext   =   set_symfuncs(@sym_ext_fhn,'sys_tau',@()[ind.tau1, ind.tau2],'sys_cond', @(p)condTestSys(p,ind));

varargin    =   {'extra_condition',true, 'matrix', 'sparse', ...
                    'newton_max_iterations', 15, 'print_residual_info',1};

default     =   {'degree',5,'intervals', 200,'corpar',[], 'step', 1e-3, 'indperiod',[],'stepcond',[]};

[options,pass_on]   =   dde_set_options(default,varargin,'pass_on');

%% Loading the bifurcation curve file
% We load the bifurcation curve and choose representative periodic solution
% which we will use to compute the master stability curve

load('./Data/BifurcationCurve.mat')
auxBifCurve =   BifCurve.branch;

auxSelectK =   [1/15, 1/16,0.0617, 1/20,1/30 ,1/50,1/100;...
                    1e-3,1e-3,1e-3,1e-4,1e-3,1e-4, 5e-4];

parK        =   arrayfun(@(x) x.parameter(ind.k), auxBifCurve.point);
perProbSol  =   cell(1,7);
auxInd      =   zeros(1,7);
for i=1:length(auxInd)
    auxOut                          =   find(abs(parK-auxSelectK(1,i))< auxSelectK(2,i));
    auxPerPoint                     =   auxBifCurve.point(auxOut(1));
    auxPerPoint.parameter(ind.k)    =   auxSelectK(1,i);

    [auxPerPoint,~]                 =   p_correc(funcs_ext,auxPerPoint,[ind.l1, ind.tau1,ind.tau2,ind.phi],options.stepcond,auxBifCurve.method.point);
    perProbSol{i}                   =   auxPerPoint;
    auxInd(i)                       =   auxOut(1);
end

%% Computing the master stability curve that contains the trivial bundle
% Here, we loop over the selected solutions in the previous section
% and compute the corresponding master stability curves

masterCurves    =   cell(1,length(auxInd));
auxBound        =   [pi,pi,pi,pi,pi,pi/50,pi/100];
auxMaxStep      =   [1.1,1.1,1.1,1.1,1.1,1/10,1/10];
for i=1:length(perProbSol)
    auxBifCurve1                    =   BifCurve.branch;
    auxBifCurve1.point(auxInd(i))   =   perProbSol{i};

    disp(1/auxBifCurve1.point(auxInd(i)).parameter(ind.k))
    [br1,suc]                       =   SetupPsolFrom_psol(funcs_ext, auxBifCurve1, ...
        auxInd(i),'contpar', [ind.l1, ind.l2, ind.phi, ind.tau1, ind.tau2], 'extra_condition','true', 'step', 1e-3);

    br1.parameter.min_bound         =   [ind.phi, -auxBound(i)];
    br1.parameter.max_bound         =   [ind.phi, auxBound(i)];
    br1.parameter.max_step          =   [0, 0.5; ind.phi, 0.01*auxMaxStep(i)];

    figure(i); clf; hold on;
    br1                             =   br_contn(funcs_ext,br1,1000);
    br1                             =   br_rvers(br1);
    br1                             =   br_contn(funcs_ext,br1,1000);
    hold off;

    masterCurves{i}.branch           =   br1;
    masterCurves{i}.k                =   br1.point(1).parameter(ind.k);
end


%% Computing secondary master stability curves if present
% Here, we compute the second stability curve that detaches from the trivial 
% one for smaller k. We perform numerical continuation to k=1/50 and k=1/100
% and compute the secondary curve

auxBifCurve1                    =   BifCurve.branch;
auxBifCurve1.point(auxInd(i))   =   perProbSol{1};
[br1,suc]                       =   SetupPsolFrom_psol(funcs_ext, auxBifCurve1, auxInd(1),'contpar', [ind.l1, ind.l2, ind.phi, ind.tau1, ind.tau2], 'extra_condition','true', 'step', 1e-3);

br1.parameter.min_bound         =   [ind.phi, -3*pi/15];
br1.parameter.max_bound         =   [ind.phi, 3*pi/15];
br1.parameter.max_step          =   [0, 1; ind.phi, 0.1];

figure(1); clf; hold on;
br1                             =   br_contn(funcs_ext,br1,1000);
hold off;

brStarting                      =   {};
[br2,suc]                       =   SetupPsolFrom_psol(funcs_ext, br1, length(br1.point),'contpar', [ind.k, ind.tau1,ind.l1, ind.l2,  ind.tau2], 'extra_condition','true', 'step', 1e-3);

br2.parameter.min_bound         =   [ind.k, 1/50];
br2.parameter.max_bound         =   [ind.k, 1/12];
br2.parameter.max_step          =   [0, 0.5; ind.k, 0.01*auxMaxStep(i)];

figure(1); clf; hold on;
br2                             =   br_rvers(br2);
br2                             =   br_contn(funcs_ext,br2,1000);
hold off;
brStarting{1}                   =   br2;

br2.parameter.min_bound         =   [ind.k, 1/100];
figure(1); clf; hold on;
br2                             =   br_contn(funcs_ext,br2,1000);
hold off;
brStarting{2}                   =   br2;

for i=1:2
    [br3,suc]                       =   SetupPsolFrom_psol(funcs_ext, brStarting{i}, ...
        length(brStarting{i}.point),'contpar', [ind.l1, ind.l2, ind.phi, ind.tau1, ind.tau2], 'extra_condition','true', 'step', 1e-3);

    br3.parameter.min_bound         =   [ind.phi, -1.1*pi];
    br3.parameter.max_bound         =   [ind.phi,  1.1*pi];
    br3.parameter.max_step          =   [0, 0.5; ind.phi, 0.005];

    figure(i); clf; hold on;
    br3                             =   br_contn(funcs_ext,br3,1000);
    br3                             =   br_rvers(br3);
    br3                             =   br_contn(funcs_ext,br3,1000);
    hold off;

    masterCurves{i+5}.branch2        =   br3;
end

%% Saving the curves in different files
auxLabels   =   {'15','16','BIF','20','30','50','100'};
for i=1:length(masterCurves)
    masterCurve         =    masterCurves{i};
    masterCurve.desc    =   'Master Stability Curve of a periodic solution along the principal branch';
    save(['./Data/MasterStabilityCurves_inv_q_' auxLabels{i}  '.mat'], 'masterCurve');
end