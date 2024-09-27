%% Script to generate the Floquet multipliers and Master Stability Curve of selected two-pulse periodic solution
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
% In this script we load the Master stability curve file, and then select 
% the periodic solution and compute their Floquet Exponent by simulating
% its monondromy matrix and finding its Eigenvalues.


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

funcs_ext   =   set_symfuncs(@sym_ext_fhn,'sys_tau',@()[ind.tau1, ind.tau2],'sys_cond', @(p)condTestSys(p,ind));

%% Computing the case 
load('./Data/BifurcationCurveTwoPulses.mat')
selectedProb    =   200;
auxM            =   1;


floqStructure.perDDE        = BifCurveTwoPulses.branch.point(end); %The last point correspond to 1/k=200
floqStructure.k             = floqStructure.perDDE.parameter(ind.k);
floqStructure.valuesM       = cell(1,1);

tic
valM                                = 1;
[floqM, floqE, auxVec]              = compFloquetExpMonBiDi(floqStructure.perDDE,valM,1e-01);
floqStructure.valuesM{1}.M          = valM;
floqStructure.valuesM{1}.floqM      = floqM;
floqStructure.valuesM{1}.floqE      = floqE;
% This is the Floquet Bundle before the exponential scaling. That is, the first point and the last point are different by factor given by the Floquet Multiplier
floqStructure.valuesM{1}.floqBun    = auxVec;
toc

save(['./Data/FloquetStructureTwoPulsePeriodic_inv_q_' num2str(selectedProb(1))  '.mat'], 'floqStructure');

%% Computing the master stability curves
load('./Data/FloquetStructureTwoPulsePeriodic_inv_q_200.mat')

masterCurveDoublePulse              =   cell(1);
masterCurveDoublePulse{1}.k         =   floqStructure.perDDE.parameter(ind.k);
masterCurveDoublePulse{1}.branches  =   cell(1);

iterLimits  =   [250,250,3200];
auxPer      =   floqStructure.perDDE;
for i=1:3
    auxM        =   floqStructure.valuesM{1}.floqM(i);
    auxVec      =   floqStructure.valuesM{1}.floqBun(:,i);
    [auxBranch] =   initBundleFromEigVector(auxPer,auxVec,auxM,ind,funcs_ext);
    close all

    [auxBranch1,suc]    =   SetupPsolFrom_psol(funcs_ext, auxBranch, length(auxBranch.point),'contpar', [ind.l1, ind.l2,ind.phi,ind.tau1, ind.tau2], 'extra_condition','true', 'step', 1e-2);
    auxBranch1.parameter.max_bound =   [ind.phi, pi]; 
    auxBranch1.parameter.min_bound =   [ind.phi, -pi];
    auxBranch1.parameter.max_step  =   [0, 1e-02 ind.phi, 0.01];
    
    figure(i); clf; hold on;
    auxBranch1  =   br_contn(funcs_ext,auxBranch1,iterLimits(i));
    auxBranch1  =   br_rvers(auxBranch1);
    auxBranch1  =   br_contn(funcs_ext,auxBranch1,iterLimits(i));
    hold off;

    masterCurveDoublePulse{1}.branches{i}   =   auxBranch1;
end
save(['./Data/MasterStabilityCurvesTwoPulse_inv_q_' num2str(selectedProb(1))  '.mat'], 'masterCurveDoublePulse');

%% Compute the Spatio-Temporal plot
load('./Data/FloquetStructureTwoPulsePeriodic_inv_q_200.mat')

auxPer = floqStructure.perDDE;

spatioTwoPulseTemporal.perDDE       =   floqStructure.perDDE;
spatioTwoPulseTemporal.k            =   spatioTwoPulseTemporal.perDDE.parameter(ind.k);
auxM                                =   1;
[solAux]                            =   compSimulationBiDi(spatioTwoPulseTemporal.perDDE,auxM,4000,0,[],[]);
spatioTwoPulseTemporal.simul        =   cell(1);
spatioTwoPulseTemporal.simul{1}.M   =   auxM;
spatioTwoPulseTemporal.simul{1}.sol =   solAux;

save(['./Data/SimulationTwoPulsePeriodic_inv_q_200.mat'], 'spatioTwoPulseTemporal');

clc
clear all