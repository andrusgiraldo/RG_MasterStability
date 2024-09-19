%% Script to generate the Spatio-Temporal plots
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
% In this script we load the Master stability curve file, and then select 
% the periodic solution and a simulate bi-directional ring network using as
% initial conditions suitable points along the periodic solution.


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

%% Computing the case for q=1/20

selectedProb        =   {15, 20, 50 ,100};
selectedMforProb    =   containers.Map(selectedProb,{[1,4],[1,5],[1,2],[1,2]});
auxShift            =   [0,0.25,0,0];
auxTf               =   [2500, 100,4000,4000];

for i=1:length(selectedProb)
    load(['./Data/MasterStabilityCurves_inv_q_' num2str(selectedProb{i})  '.mat'])
    spatioTemporal.perDDE        =  masterCurve.branch.point(1);
    spatioTemporal.k             =  spatioTemporal.perDDE.parameter(ind.k);
    auxM                         =  selectedMforProb(selectedProb{i});
    for j=1:length(auxM)
        tic
        auxEmtpy =[];
        if(selectedProb{i}==15 && auxM(j)==4)
            % Perturbation in scale for Figure 2 panel (d)
            rng(4);     pert = 1 + 0.3*rand(1,15*4);
            [solAux]                      = compSimulationBiDi(spatioTemporal.perDDE,auxM(j),auxTf(i),auxShift(i),pert,auxEmtpy);
        elseif(selectedProb{i}==15 && auxM(j)==1)
            rng(1);     pert = 0.7 + 0.3*rand(1,15*4);
            [solAux]                      = compSimulationBiDi(spatioTemporal.perDDE,auxM(j),auxTf(i),auxShift(i),pert,auxEmtpy);
        else
            [solAux]                      = compSimulationBiDi(spatioTemporal.perDDE,auxM(j),auxTf(i),auxShift(i),[],auxEmtpy);
        end
        spatioTemporal.simul{j}.M     = auxM(j);
        spatioTemporal.simul{j}.sol   = solAux;
        toc
    end
    save(['./Data/SimulationMainPeriodic_inv_q_' num2str(selectedProb{i})  '.mat'], 'spatioTemporal');
end