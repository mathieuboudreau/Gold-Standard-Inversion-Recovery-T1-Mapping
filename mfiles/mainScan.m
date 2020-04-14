% mainScan.m
%
% Loads .mat file from the directory 'data' in the T1path, 
% performs T1 mapping, and displays the results 
%
% written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%  (c) Board of Trustees, Leland Stanford Junior University

clear all
close all

T1path = '../';

%% Where to find the data
loadpath = [T1path 'data/'];

datasetnb = 2;

switch (datasetnb)
	case 1
		filename = 'TestSingleSlice'; % complex fit
		method = 'RD-NLS'
	case 2
		filename = 'IRData.mat'; % magnitude fit
		method = 'RD-NLS-PR'
end

%% Properties of the data
%

extra.T1Vec = 1:5000;  % Initial grid points for the T1 search
extra.tVec = [350, 500, 650, 800, 950, 1100, 1250, 1400, 1700]; % Inversion times (TIs) considered

%% Where to save the data
savepath = [T1path 'fitdata/'] 

loadStr = [loadpath filename]
saveStr = [savepath 'T1Fit' method '_' filename]

%% Perform fit
T1ScanExperiment(loadStr, saveStr, method, extra);

%% Display results
T1FitDisplayScan(loadStr, saveStr);
