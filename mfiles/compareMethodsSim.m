% compareMethodsSim.m
%
% Performs a Monte-Carlo simulation of T1 mapping experiments. Different
% methods are compared on the same dataset. Specific parameters are 
% specified below.
%
% This script is used to produce Table 1 in the article.
%
% written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%  (c) Board of Trustees, Leland Stanford Junior University

clear all
close all

saveData = 1; % 1 - save data, 0 - no save
dofit = 1; % 1 - do the fit, 0 = reuse saved simulation data
dodisplay = 1; % 1 - display results, 0 - no display

%% Set the algorithms to be compared
typeOfSim = {'LM','LMSamePhase','RD-NLS',...
  'LMSamePhaseMag','LMSamePhasePR','RD-NLS-PR'}
%typeOfSim = {'RD-NLS','RD-NLS-PR'}

nbrOfSims = length(typeOfSim);
%% Parameters 
T1 = 263       % True T1
MC = 20000       % Number of Monte-Carlo simulations, 20000
flipAngle = 172;  % Effective flip angle

extra.TR = 2550;  % Repetition time (TR)
extra.T1Vec = 1:5000;  % Initial grid points for the T1 search
extra.tVec = [50,400,1100,2500]; % Inversion times (TIs) considered 
extra.theta = flipAngle;

M0 = 1; % Normalized signal
stdNoise = 0.03;

%% Generate the data
tLen = length(extra.tVec);

% Noise-free signal
S_noNoise = M0*( 1 - (1-cos(extra.theta*pi/180))*exp(-extra.tVec/T1) - cos(extra.theta*pi/180)*exp(-extra.TR/T1) ).';
% Add noise
noise = stdNoise/sqrt(2).*(randn(tLen,MC)+i*randn(tLen,MC));
S = repmat(S_noNoise,1,MC) + noise;

%% Fit %%

if (dofit)
  disp(['Processing ' num2str(MC) ' voxels, using ' ...
    num2str(nbrOfSims) ' types of algorithms.'])
  startTimeTotal = cputime;
  timeTaken = zeros(nbrOfSims,1);
  T1Est = zeros(MC,nbrOfSims);
  res = zeros(MC,nbrOfSims);
  for k = 1:nbrOfSims
    switch(typeOfSim{k})
      case 'LM'
        % Initial value for LM fit
        extra.kInit =  2; extra.T1Init =  500; 
        startTimeMC = cputime;
        for mc = 1:MC
          % Display progress on screen
          if rem(mc,1000) == 0
            disp(['LM, at mc=' num2str(mc)])
          end
          [T1Est(mc,k)] = lm(S(:,mc), extra);
        end
        timeTaken(k) = cputime - startTimeMC;
        disp(['Processed ' typeOfSim{k} ' in ' ...
          num2str(timeTaken(k)) ' sec.'])
      case 'LMSamePhase'
        % Initial value for LM same phase fit
        extra.kInit =  2; extra.T1Init =  500;
        startTimeMC = cputime;
        for mc = 1:MC
          [T1Est(mc,k)] = lmSph(S(:,mc), extra);
        end
        timeTaken(k) = cputime - startTimeMC;
        disp(['Processed ' typeOfSim{k} ' in ' ...
          num2str(timeTaken(k)) ' sec.'])
      case 'RD-NLS'
        nlsS = getNLSStruct(extra);
        startTimeMC = cputime;
        for mc = 1:MC
          [T1Est(mc,k)] = rdNls(S(:,mc), nlsS);
        end
        timeTaken(k) = cputime - startTimeMC;
        disp(['Processed ' typeOfSim{k} ' in ' ...
          num2str(timeTaken(k)) ' sec.'])
      case 'LMSamePhaseMag'
        % Initial value for LM same phase fit with magnitude data
        extra.kInit =  2; extra.T1Init =  500;
        startTimeMC = cputime;
        for mc = 1:MC
          % Display progress on screen
          if rem(mc,1000) == 0
            disp(['Mag LM, at mc=' num2str(mc)])
          end
          [T1Est(mc,k)] = lmSphMag(abs(S(:,mc)).^2, extra);
        end
        timeTaken(k) = cputime - startTimeMC;
        disp(['Processed ' typeOfSim{k} ' in ' ...
          num2str(timeTaken(k)) ' sec.'])
      case 'LMSamePhasePR'
        % Initial value for LM same phase fit with magnitude data and
        % polarity restoration
        extra.kInit =  2; extra.T1Init =  500;
        startTimeMC = cputime;
        for mc = 1:MC
          % Display progress on screen
          if rem(mc,1000) == 0
            disp(['Mag LM PR, at mc=' num2str(mc)])
          end
          [T1Est(mc,k)] = lmSphPr(abs(S(:,mc)), extra);
        end
        timeTaken(k) = cputime - startTimeMC;
        disp(['Processed ' typeOfSim{k} ' in ' ...
          num2str(timeTaken(k)) ' sec.'])
      case 'RD-NLS-PR'
        nlsS = getNLSStruct(extra);
        startTimeMC = cputime;
        for mc = 1:MC
          [T1Est(mc,k)] = rdNlsPr(abs(S(:,mc)), nlsS);
        end
        timeTaken(k) = cputime - startTimeMC;
        disp(['Processed ' typeOfSim{k} ' in ' ...
          num2str(timeTaken(k)) ' sec.'])
    end
  end
  timeTakenTot = cputime - startTimeTotal;
  disp(['Total time taken: ' num2str(timeTakenTot/60) 'min'])
else
  load(['../fitdata/simDataCompareAllMethods'  ...
    '_MC' num2str(MC) '_' num2str(T1) 'T1'])
  dodisplay = 1;
end
clear S mc k noise startTimeMC startTimeTotal
clear aEst bEst aMagEst bMagEst kEst cEst cMagEst

%% Display the results
if (dodisplay)
  % Plot the histograms
  for k = 1:nbrOfSims
    titleTmp = [typeOfSim{k} ', T1 for \sigma_n=' num2str(stdNoise)];
    T1FitDisplaySim(round(T1Est(:,k)),extra,titleTmp);
  end
  clear k myMu mySigma tmp1 tmp2 tmp newTitle
end

%% Save the results
if (saveData)
  clear saveData
  save(['../fitdata/simDataCompareAllMethods_MC' num2str(MC) '_' num2str(T1) 'T1'])
end
