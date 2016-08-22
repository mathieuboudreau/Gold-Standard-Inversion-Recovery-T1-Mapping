% [T1Est, bEst, aEst, res] = lm(data, extra)
%
% Finds estimates of T1, a, and b using the Levenberg-Marquardt
% algorithm via fminsearch. The model a+b*exp(-t/T1) is used. 
% The residual is the rms error between the data and the fit. 
% 
% INPUT:
% data - the data to estimate from
% extra.tVec - vector of TI's used
%
% The initialization of the LM algorithm can be done in two ways:
% 1)
% extra.x0   - 5x1 vector of starting points, such that
%              a  = x(1) + i*x(2)
%              b  = x(3) + i*x(4)
%              T1 = x(5).
% OR 
% 2)
% extra.T1Init 
% extra.kInit (typically set to 2 as k = 1-cos(flipangle)). 
% The other unknowns are given as the the data at the largest TI.
%    
% written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%  (c) Board of Trustees, Leland Stanford Junior University 

function [T1Est, bEst, aEst, res] = lm(data, extra)
    
%% Set the initial values for the search
[tmp,order] = sort(extra.tVec);
try
  x0(5) = extra.T1Init;
  x0(1) = real(data(order(end)));
  x0(2) = imag(data(order(end)));
  x0(3) = -extra.kInit*x0(1);
  x0(4) = -extra.kInit*x0(2);
catch
  x0 = extra.x0;
end

%% Make sure data is a column vector
data = data(:);
extra.tVec = extra.tVec(:);

%% Do the fit
x = fminsearch( ...
  @(x)sum(abs( data-( (x(1)+i*x(2)) + (x(3)+i*x(4))*exp(-extra.tVec/x(5)) ) ).^2), ...
  x0,optimset('display','off'));

aEst = x(1) + i*x(2);
bEst = x(3) + i*x(4);
T1Est = x(5);

%% Compute the residual
modelValue = aEst + bEst*exp(-extra.tVec/T1Est); 
res = 1/sqrt(length(data))*norm(1 - modelValue./data);
