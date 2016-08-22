% [T1Est, kEst, cEst, res] = lmSphPr(data, extra)
%
% Finds estimates of T1, |c|, and k=-1+cos(theta) using the 
% Levenberg-Marquardt algorithm via fminsearch and phase restoration.
% The model |c|*(1-k*exp(-t/T1)) is used, i.e., there is only one phase
% and we have taken the magnitude of the data.
% The residual is the rms error between the data and the fit. 
% 
% INPUT:
% data - the data to estimate from
% extra.tVec - vector of TI's used
%
% The initialization of the LM algorithm can be done in two ways:
% 1)
% extra.x0   - 3x1 vector of starting points, such that
%              |c| = x(1)
%              k  = x(2)
%              T1 = x(3).
% OR 
% 2)
% extra.T1Init 
% extra.kInit (typically set to 2 as k = 1-cos(flipangle)). 
% |c| is given as the the data at the largest TI.
%    
% written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%  (c) Board of Trustees, Leland Stanford Junior University 

function [T1Est, kEst, cMagEst, res] = lmSphPr(data, extra)
   
%% Set the initial values for the search
% Make sure the data come in increasing TI-order
[extra.tVec,order] = sort(extra.tVec); 
data = squeeze(data); 
data = data(order);

try
  x0(3) = extra.T1Init;
  x0(2) = extra.kInit;
  x0(1) = abs(data((end)));
catch
  x0 = extra.x0;
end

%% Make sure data vectors are a column vectors
data = data(:);
extra.tVec = extra.tVec(:);
N = length(extra.tVec);
x = zeros(3,2);
resTmp = zeros(1,2);

%% Find the min of the data
[minVal, minInd] = min(data);

%% Fit
for ii = 1:2
  if ii == 1
    % First, we set all elements up to and including
    % the smallest element to minus
    dataTmp = data.*[-ones(minInd,1); ones(N - minInd,1)];
  elseif ii == 2
    % Second, we set all elements up to (not including)
    % the smallest element to minus
    dataTmp = data.*[-ones(minInd-1,1); ones(N - (minInd-1),1)];
  end
  
  % Do the fit
  x(:,ii) = fminsearch( ...
    @(x)sum( ( dataTmp - ...
    x(1)*(1 - x(2)*exp(-extra.tVec/x(3))) ).^2 ), ...
    x0,optimset('display','off'));
  resTmp(ii) = norm( ...
    1 - x(1,ii)*(1 - x(2,ii)*exp(-extra.tVec/x(3,ii)))./dataTmp);
end
  
[res,ind] = min(resTmp);
res = res/sqrt(N);
cMagEst = x(1,ind);
kEst = x(2,ind);
T1Est = x(3,ind);
