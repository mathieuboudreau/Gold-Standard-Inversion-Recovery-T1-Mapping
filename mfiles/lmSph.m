% [T1Est, kEst, cEst, res] = lmSph(data, extra)
%
% Finds estimates of T1, c, and k=-1+cos(flipangle) using the 
% Levenberg-Marquardt algorithm via fminsearch. 
% The model c(1-k*exp(-t/T1)) is used, i.e. there is only one phase. 
% The residual is the rms error between the data and the fit. 
% 
% INPUT:
% data - the data to estimate from
% extra.tVec - vector of TI's used
%
% The initialization of the LM algorithm can be done in two ways:
% 1)
% extra.x0   - 4x1 vector of starting points, such that
%              c  = x(1) + i*x(2)
%              k  = x(3)
%              T1 = x(4).
% OR 
% 2)
% extra.T1Init 
% extra.kInit (typically set to 2 as k = 1-cos(flipangle)). 
% c is given as the the data at the largest TI.
%    
% written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%  (c) Board of Trustees, Leland Stanford Junior University 

function [T1Est, kEst, cEst, res] = lmSph(data, extra)
    
%% Set the initial values for the search
[tmp,order] = sort(extra.tVec);
try
  x0(4) = extra.T1Init;
  x0(3) = extra.kInit;
  x0(1) = real(data(order(end)));
  x0(2) = imag(data(order(end)));
catch
  x0 = extra.x0;
end

%% Make sure data is a column vector
data = data(:);
extra.tVec = extra.tVec(:);

%% Do the fit
x = fminsearch( ...
  @(x)sum(abs( data-( (x(1)+i*x(2))*(1 - x(3)*exp(-extra.tVec/x(4))) ) ).^2), ...
  x0,optimset('display','off'));

cEst = x(1) + i*x(2);
kEst = x(3);
T1Est = x(4);

%% Compute the residual
modelValue = cEst*(1 - kEst*exp(-extra.tVec/T1Est)); 
res = ...
  1/sqrt(length(data))*norm(1 - modelValue./data);
