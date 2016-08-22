% mainSim.m
%
% Performs a Monte-Carlo simulation of a T1 mapping experiment. 
% Specific parameters are specified below. 
%
% written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%  (c) Board of Trustees, Leland Stanford Junior University

clear all
close all


%Reset randn to initial state...
randn('state',0)
%...or make sure it's reset every time
%randn('state',sum(100*clock))

saveData = 1; %1 - save data, 0 - no save
dofit = 1; % 1 - do the fit, 0 = reuse saved simulation data
dodisplay = 1; % 1 - display results, 0 - no display

%% Set the method to test
method = 'RD-NLS' % Using the fast RD-NLS search for complex data
%method = 'RD-NLS-PR' % Fast search for magnitude data 

%method = 'LM' % Using fminsearch with five parameters
%method = 'LMSamePhase' % Using fminsearch with model having same phase

%method = 'LMSamePhaseMag' % Using fminsearch with model having same phase and with magnitude-squared data
%method = 'LMSamePhasePR' % Using fminsearch with model having same phase and with magnitude data using a phase restoration method

%% Parameters 
T1 = 263;         % True T1
stdNoise = 0.03; 
MC = 20000;        % Number of Monte-Carlo simulations %20000
flipAngle = 172;  % Effective flip angle
extra.TR = 2550;  % Repetition time (TR)
extra.T1Vec = 1:5000;  % Initial grid points for the T1 search
extra.tVec = [50,400,1100,2500]; % Inversion times (TIs) considered

%% Fit 
if (dofit)
	% The LM methods require initial values. They are set here.
	switch(method)
		case{'LM'}
			%extra.x0 =  [10 10 2 2 200];
			extra.kInit =  2;
			extra.T1Init =  200;
		case{'LMSamePhase'}
			%Use this initialization...
			%extra.x0 = [1 1 2 200];
			
			%...or this
			extra.kInit =  2; extra.T1Init =  200;
		case{'LMSamePhaseMag'}
			%Use this initialization...
			%extra.x0 = [1 2 200];
			
			%...or this
			extra.kInit =  2; extra.T1Init =  200;
		case{'LMSamePhasePR'}
			%Use this initialization...
			%extra.x0 = [1 2 200];
			
			%...or this
			extra.kInit =  2; extra.T1Init =  200;
	end
	[T1Est, bEst, aEst, res] = ...
		T1SimExperiment(MC,stdNoise,T1,flipAngle,extra,method);
else
	load(['../fitdata/T1Sim_' method '_MC' num2str(MC) '_' num2str(T1) 'T1']);
	dodisplay = 1;
end
  
%% Display results
if (dodisplay)
	histTitle = ['T1, ' method ', \sigma_n=' num2str(stdNoise) ];
	T1FitDisplaySim(T1Est,extra,histTitle);
end

%% Save results
if (saveData)
  clear saveData
  save(['../fitdata/T1Sim_' method '_MC' num2str(MC) '_' num2str(T1) 'T1'])
end

