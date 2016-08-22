% T1ScanExperiment(loadStr, saveStr, method)
%
% loadStr: the data on which the fit will be performed, to be loaded
% saveStr: where results of the fit will be saved
% method:  what fitting method to use: RD-NLS, RD-NLS-PR
%
% Estimates T1 together with:
%   RD-NLS: a and b parameters to fit the data to a + b*exp(-TI/T1)
%   RD-NLS-PR: ra and rb parameters to fit the data to |ra + rb*exp(-TI/T1)|                                                  
%
% written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%  (c) Board of Trustees, Leland Stanford Junior University

function T1ScanExperiment(...
  loadStr, saveStr, method)

% Load the data on which to perform the fit
load(loadStr);

savefitdata = 1; 

nlsS = getNLSStruct(extra,1);

switch method
  case{'RD-NLS'}
    nbrOfFitParams = 4; % Number of output arguments for the fit
    fitStr = ['[T1Est bEst aEst res] = rdNls(data(jj,:),nlsS);']
    storeStr = 'll_T1(jj,:) = [T1Est bEst aEst res];'
    clearStr = 'clear jj T1Est bEst aEst res'
  case{'RD-NLS-PR'}
    data = abs(data);
    nbrOfFitParams = 4; % Number of output arguments for the fit
    fitStr = ['[T1Est bMagEst aMagEst res] = rdNlsPr(data(jj,:),nlsS);']
    storeStr = 'll_T1(jj,:) = [T1Est bMagEst aMagEst res];'
    clearStr = 'clear jj T1Est bMagEst aMagEst res'
end

dims = size(data);
nbrow = size(data,1);
nbcol = size(data,2);

if numel(dims) > 3
    nbslice = dims(3); % Check number of slices 
else
    nbslice = 1;
    tmpData(:,:,1,:) = data; % Make data a 4-D array regardless of number of slices
    data = tmpData;
    clear tmpData;
end

% The data needs to lie along an exponential.
% Refer to getData.m if that is not the case.
while (true)
	if (nbslice ~= 1)
		disp('For which slice would you like to check the data?');
		zz = input(['Enter number 1 to ' num2str(nbslice) '.  0 for no check --- '], 's');
		zz = cast(str2num(zz), 'int16');
		if (isinteger(zz) & zz >= 0 & zz <= nbslice) break;
		end
	else
		zz = input(['Enter 1 to check the data, 0 for no check --- '], 's');
		zz = cast(str2num(zz), 'int16');
		if (isinteger(zz) & zz >= 0 & zz <= nbslice) break;
		end
	end	
end

if (zz ~= 0)
	sliceData = squeeze(data(:,:,zz,:));
	disp('Click on one point to check that the time series looks like an exponential')
	disp('Refer to getData.m if that is not the case.  CTRL-click or right-click when done');
	TI = extra.tVec;
	plotData(real(sliceData),TI);
end

close all

dataOriginal = data;

% Mask the background (points dimmer than maskFactor*the brightest point) 
% Done on the last image, where we expect the magnetization to have almost
% fully recovered. 
% Increase mF is the mask does not cover all the background. Decrease it if
% it covers part of the object. 
mF = 0.1;
maskFactor = mF;
mask = zeros(nbrow, nbcol, nbslice);

[u,v] = max(extra.tVec);

for kk = 1:nbslice
	maskTmp = mask(:,:,kk);
	maskTmp = medfilt2(maskTmp); % remove salt and pepper noise
	maskThreshold = maskFactor*max(max(abs(data(:,:,kk,v))));
	maskTmp(find(abs(data(:,:,kk,v))> maskThreshold)) = 1;
	mask(:,:,kk) = maskTmp;
	clear maskTmp
end

figure, 
imshow(abs(mask(:,:,1)),[])

disp('Check that maskFactor is appropriate. Hit enter if so. Cancel and adjust mF in T1ScanExperiment.m if not.')
pause

maskInds = find(mask);

nVoxAll = length(maskInds);
% How many voxels to process before printing out status data
numVoxelsPerUpdate = min(floor(nVoxAll/10),1000); 
						   
ll_T1 = zeros(nVoxAll, nbrOfFitParams);
% Number of status reports
nSteps = ceil(nVoxAll/numVoxelsPerUpdate); 

for ii = 1:size(data,4)
    tmpVol = data(:,:,:,ii);
    tmpData(:, ii) = tmpVol(maskInds)';
end
clear tmpVol;
data = tmpData;
clear tmpData;

startTime = cputime;
fprintf('Processing %d voxels.\n', nVoxAll);

h = waitbar(0, sprintf('Processing %d voxels', nVoxAll)); 

for ii = 1:nSteps
  curInd = (ii-1)*numVoxelsPerUpdate+1;
  endInd = min(curInd+numVoxelsPerUpdate,nVoxAll);
  for jj = curInd:endInd
    % Do the fit
    eval(fitStr) 
    % Store the data
    eval(storeStr);
  end
  waitbar(ii/nSteps, h, sprintf('Processing %d voxels, %g percent done...\n',nVoxAll,round(endInd/nVoxAll*100)));
end
eval(clearStr)
close(h);
timeTaken = round(cputime - startTime);
fprintf('Processed %d voxels in %g seconds.\n',nVoxAll, timeTaken);

dims = [size(mask) 4];
im = zeros(size(mask));

for ii = 1:nbrOfFitParams
    im(maskInds) = ll_T1(:,ii);
    T1(:,:,:,ii) = im;
end

% Going back from a numVoxels x 4 array to nbrow x nbcol x nbslice
ll_T1 = T1;

% Store ll_T1 and mask in saveStr
% For the complex data, ll_T1 has four parameters 
% for each voxel, namely:
% (1) T1 
% (2) 'b' or 'rb' parameter 
% (3) 'a' or 'ra' parameter
% (4) residual from the fit

if (savefitdata)
	save(saveStr,'ll_T1','mask','nlsS')
end

% Check the fit
TI = extra.tVec;
nbtp = 20;
timef = linspace(min(TI),max(TI),nbtp);

% Inserting a short pause, otherwise some computers seem
% to get problems
pause(1);

zz = 0;
while(true)
  if (nbslice ~= 1)
    disp('For which slice would you like to check the fit?');
    zz = input(['Enter number 1 to ' num2str(nbslice) '.  0 for no check --- '], 's');
    zz = cast(str2num(zz), 'int16');
    if (isinteger(zz) & zz >= 0 & zz <= nbslice) 
      break;
    end
  else
    zz = input(['Enter 1 to check the fit, 0 for no check --- '], 's');
    zz = cast(str2num(zz), 'int16');
    if (isinteger(zz) & zz >= 0 & zz <= nbslice) 
      break;
    end
  end
end

if (zz ~= 0)
  sliceData = squeeze(dataOriginal(:,:,zz,:));
  datafit = zeros(nbrow,nbcol,nbtp);
  switch method
    case{'RD-NLS'}
		for kk = 1:nbtp
			datafit(:,:,kk) = ll_T1(:,:,zz,3) + ...
				ll_T1(:,:,zz,2).*exp(-timef(kk)./ll_T1(:,:,zz,1));
		end
	case{'RD-NLS-PR'}
		for kk = 1:nbtp
			datafit(:,:,kk) = abs(ll_T1(:,:,zz,3) + ...
				ll_T1(:,:,zz,2).*exp(-timef(kk)./ll_T1(:,:,zz,1)));
		end
  end
  disp('Click on one point to check the fit. CTRL-click or right-click when done')
  plotData(real(sliceData),TI,real(datafit),ll_T1);
  close all
end
