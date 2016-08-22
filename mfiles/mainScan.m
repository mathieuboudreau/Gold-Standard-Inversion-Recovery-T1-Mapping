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

datasetnb = 1;

switch (datasetnb)
	case 1
		filename = 'TestSingleSlice'; % complex fit
		method = 'RD-NLS'
	case 2
		filename = 'TestSingleSlice'; % magnitude fit
		method = 'RD-NLS-PR'
end


%% Where to save the data
savepath = [T1path 'fitdata/'] 

loadStr = [loadpath filename]
saveStr = [savepath 'T1Fit' method '_' filename]

%% Perform fit
T1ScanExperiment(loadStr, saveStr, method);

%% Display results
T1FitDisplayScan(loadStr, saveStr);
