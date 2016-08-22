% [sigma,mu] = T1FitDisplaySim(T1Est,extra,myTitle)
%
% written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%  (c) Board of Trustees, Leland Stanford Junior University

function [sigma,mu] = T1FitDisplaySim(T1Est, extra, myTitle)

T1EstTmp = T1Est(:);
T1EstTmp = T1EstTmp(find(T1EstTmp > extra.T1Vec(1)));
T1EstTmp = T1EstTmp(find(T1EstTmp < extra.T1Vec(end)));
T1Est = T1EstTmp;

x = min(T1Est(:)):1:max(T1Est(:));
x = x';
h = histc(T1Est(find(T1Est>0)),x); % centered on integers
%hs = medfilt1(h,10); % when low SNR, 10ms median filter cf. Tozer03 or Tofts 

hs = h; 
[sigma,mu,A] = customGaussFit(x,hs);

y = A*exp(-(x-mu).^2/(2*sigma^2));

figure
bar(x,hs,1)
hold
plot(x,y,'.r')
customFormat

if (nargin == 3)
  if isempty(myTitle)
    title([' Peak at ',num2str(round(mu)),...
      ' ms, \sigma = ' num2str(round(sigma))])
  else
    title([myTitle, 10,...
      ' Peak at ',num2str(round(mu)),...
      ' ms, \sigma = ' num2str(round(sigma))])
  end
else
  title(['T1', 10,...
    ' Peak at ',num2str(round(mu)),...
    ' ms, \sigma = ' num2str(round(sigma))])
end
  