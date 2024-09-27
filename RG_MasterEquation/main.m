%% Main execution file to produce data and figures in [Ref]
% By Andrus Giraldo and Stefan Ruschel - 03/09/2024

clc
clear all

% Generate the bifurcation curve for the main pulse
genBifCurve

% Generate the master stability curve for selected periodic profiles
genMasterCurves

% Compute the Floquet exponents (also multipliers) of selected periodic
% problems by computing the Monondromy matrix
genFloquetPoints

% Create selected Spatio Temporal plots
genSTPData

% Generate the bifurcation curve for the two pulse
genBifCurveTwoPulses

% Generate the master stability curve for selected periodic profile, also
% its Floquet exponents by computing its monondromy matrix. Also, the
% spatio temporal plot of a selected profile is computed
genFloquetMasterTwoPulses

% Create all the figures of the paper
mainFigGen