% [sigma,mu,A] = customGaussFit(x,y,h)
%
% this function is doing fit to the function
% y = A * exp( -(x-mu)^2 / (2*sigma^2) )
%
% the fitting is been done by a polyfit
% the lan of the data.
%
% h is the threshold which is the fraction
% from the maximum y height that the data
% is been taken from.
% h should be a number between 0-1.
% if h have not been taken it is set to be 0.2
% as default.
%
% http://www.mathworks.com/matlabcentral/fileexchange/11733-gaussian-curve-fit
% Copyright (c) 2007, Yohanan Sivan
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

function [sigma,mu,A] = customGaussFit(x,y,h)

%% threshold
if nargin == 2, h=0.2; end

%% cutting
ymax = max(y);
xnew = [];
ynew = [];

for n = 1:length(x)
    if y(n) > ymax*h;
        xnew = [xnew,x(n)];
        ynew = [ynew,y(n)];
    end
end

%% fitting
ylog = log(ynew);
xlog = xnew;
p = polyfit(xlog,ylog,2);
A2 = p(1);
A1 = p(2);
A0 = p(3);
sigma = sqrt(-1/(2*A2));
mu = A1*sigma^2;
A = exp(A0+mu^2/(2*sigma^2));

