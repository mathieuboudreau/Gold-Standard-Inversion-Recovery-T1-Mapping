% [T1Est, bEst, aEst, res] = ...
%  T1SimExperiment(MC,stdNoise,T1,theta,extra,method)
%
% Estimates T1 for simulated experiment.
%
% written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%  (c) Board of Trustees, Leland Stanford Junior University

function [T1Est, bEst, aEst, res] = ...
  T1SimExperiment(MC, stdNoise, T1, theta, extra, method)

startTime = cputime; 

numVoxelsPerUpdate = 1000;

tLen = length(extra.tVec);
M0 = 1; % Normalized signal
extra.theta = theta;

% Noise-free signal
S_noNoise = M0*( 1 - (1-cos(theta*pi/180))*exp(-extra.tVec/T1) - ...
  cos(theta*pi/180)*exp(-extra.TR/T1) ).';

nlsS = getNLSStruct(extra);
nSteps = ceil(MC/numVoxelsPerUpdate);

fprintf('Processing %d voxels.\n', MC);
h = waitbar(0, sprintf('Simulating %d voxels', MC));

T1Est = zeros(1,MC);
bEst = zeros(1,MC);
aEst = zeros(1,MC);
res = zeros(1,MC);

for ii = 1:nSteps
  curInd = (ii-1)*numVoxelsPerUpdate+1;
  endInd = min(curInd+numVoxelsPerUpdate,MC);
  for mc = curInd:endInd
    
    %Add noise to data
    S = S_noNoise + ...
      stdNoise/sqrt(2).*(randn(tLen,1)+i*randn(tLen,1));
    
    % Do the fit
    switch(method)
      case 'RD-NLS'
        [T1Est(mc), bEst(mc), aEst(mc), res(mc)] = ...
          rdNls(S, nlsS);
      case 'RD-NLS-PR' 
        [T1Est(mc), bEst(mc), aEst(mc), res(mc)] = ...
          rdNlsPr(abs(S), nlsS);
      case 'LM'
        [T1Est(mc), bEst(mc), aEst(mc), res(mc)] = ...
          lm(S, extra);
      case 'LMSamePhase'
        [T1Est(mc), bEst(mc), aEst(mc), res(mc)] = ...
          lmSph(S, extra);
      case 'LMSamePhaseMag' 
        [T1Est(mc), bEst(mc), aEst(mc), res(mc)] = ...
          lmSphMag(abs(S).^2, extra); 
      case 'LMSamePhasePR' 
        [T1Est(mc), bEst(mc), aEst(mc), res(mc)] = ...
          lmSphPr(abs(S), extra);
    end
  end
  waitbar(ii/nSteps, h, ...
    sprintf('Processing %d voxels, %g percent done...\n',...
    MC,round(endInd/MC*100)));
end
close(h)
timeTaken = round(cputime - startTime);
fprintf('Processed %d voxels in %g seconds.\n',...
  MC, timeTaken);

