%% Script to generate the Floquet multipliers of selected periodic solution
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

%% Computing the case for q=1/20

selectedProb    =   [15, 20,50,100];
auxM            =   {[1,4],[1,5],[1,2],[1,2]};
for i=1:length(selectedProb)
    load(['./Data/MasterStabilityCurves_inv_q_' num2str(selectedProb(i))  '.mat'])
    floqStructure.perDDE        = masterCurve.branch.point(1);
    floqStructure.k             = floqStructure.perDDE.parameter(ind.k);
    floqStructure.valuesM       = cell(1,2);
    for j=1:length(auxM{i})
        tic
        valM                                = auxM{i};
        [floqM, floqE, auxVec]              = compFloquetExpMonBiDi(floqStructure.perDDE,valM(j),1);
        floqStructure.valuesM{j}.M          = valM(j);
        floqStructure.valuesM{j}.floqM      = floqM;
        floqStructure.valuesM{j}.floqE      = floqE;
        % This is the Floquet Bundle before the exponential scaling. That is, the first point and the last point are different by factor given by the Floquet Multiplier
        floqStructure.valuesM{j}.floqBun    = auxVec;
        toc
    end
    save(['./Data/FloquetStructureMainPeriodic_inv_q_' num2str(selectedProb(i))  '.mat'], 'floqStructure');
end