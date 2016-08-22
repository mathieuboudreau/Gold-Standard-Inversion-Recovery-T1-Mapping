% [T1Est, kEst, cEst, res] = lmSphMag(data, extra)
%
% Finds estimates of T1, |c|, and k=-1+cos(theta) using the 
% Levenberg-Marquardt algorithm via fminsearch. 
% The model |c*(1-k*exp(-t/T1))|^2 is used, i.e., there is only one phase
% and we have taken the magnitude-square of the data.
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
% |c| is given as the the sqrt(data) at the largest TI.
%    
% written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%  (c) Board of Trustees, Leland Stanford Junior University 

function [T1Est, kEst, cMagEst, res] = lmSphMag(data, extra)
    
%% Set the initial values for the search
[tmp,order] = sort(extra.tVec);
try
  x0(3) = extra.T1Init;
  x0(2) = extra.kInit;
  x0(1) = sqrt( abs(data(order(end))) );
catch
  x0 = extra.x0;
end

%% Make sure data is a column vector
data = data(:);
extra.tVec = extra.tVec(:);

%% Do the fit
x = fminsearch( ...
  @(x)sum( ( data - abs( x(1)*(1 - x(2)*exp(-extra.tVec/x(3))) ).^2 ).^2 ), ...
  x0,optimset('display','off'));

cMagEst = x(1);
kEst = x(2);
T1Est = x(3);

%% Compute the residual
modelValue = abs(cMagEst*(1 - kEst*exp(-extra.tVec/T1Est))).^2; 
res = ...
  1/sqrt(length(data))*norm(1 - modelValue./data);
