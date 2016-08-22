% getData.m
%
% Generates the .mat from the DICOM images for GS T1 mapping. 
%
% written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%  (c) Board of Trustees, Leland Stanford Junior University    

function getDataMag(T1path, indeces);

% clear all
% close all
% 
% T1path = '/data/mril/mril11/stikov/Scans/christine_ir_t1';

% Where to find the DICOM data 
%loadpath = [T1path 'singleslicedicomfiles/'];

loadpath = T1path;

% Where to save the .mat 
%savename = [T1path 'data/' 'TestSingleSlice'];

savename = fullfile (T1path, 'data.mat');

currpath = pwd;
cd(loadpath)
%addpath([T1path '/mfiles']);

d = mincLoadAllSeries('.');
% sized = size(d)
d = d(indeces); % Selecting only the right indexes

cd(currpath)

nbrow = size(d(1).imData,1);
nbcol = size(d(1).imData,2);
nbslice = size(d(1).imData, 3);
nbseries = length(d)

data = zeros(nbcol,nbrow,nbslice,nbseries); % Complex data
extra.tVec = zeros(1,nbseries); % One series corresponds to one TI

% % ym:
% disp('plot d(7).imData')
% figure;
% imshow(d(1).imData,[]);

for k = 1:nbseries
	dataTmp = d(k).imData;
   % sizedataTmp = size(dataTmp)
	dataTmp = double(squeeze(dataTmp));	
   % sizeSqueezeDataTmp = size(dataTmp)
	for ss = 1:nbslice
        % ym:
        % sizedata = size(data)
		%data(:,:,ss,k) = dataTmp(:,:,3+(ss-1)*4)+i*dataTmp(:,:,4+(ss-1)*4); %SE
		%data(:,:,ss,k) = dataTmp(:,:,1+(ss-1)*4); %GRE
		data(:,:,ss,k) = reshape(dataTmp(:,:,ss), nbcol, nbrow); 
        sizedatareshape = size(data)
	end
	extra.tVec(k) = d(k).inversionTime;
end 

% This is where you would correct the phase of certain datapoints. 
% Indeed, the chopping procedure during Prescan can give an additional 180?
% phase on the whole image for certain TIs (make also sure you use manual
% Prescan and don't change gains/center frequency from TI to TI!)
% e.g., data(:,:,3) = -data(:,:,3); 
%    or data(:,:,:,3) = -data(:,:,:,3); if you have more than one slice
% !!! :,:,1 is the first series acquired, not the smallest TI !!!

extra.T1Vec = 1:5000; % this range can be reduced if a priori information is available

%data = abs(data);

TI = extra.tVec
save(savename,'data','extra')

%ym:
disp('display data for the first TI');
figure;
imshow(data(:,:,1,1),[]);
hold;

disp('display data for the second TI');
figure;
imshow(data(:,:,1,2),[]);

hold;

disp('display data for the third TI');
figure;
imshow(data(:,:,1,3),[]);

hold;

disp('display data for the fourth TI');
figure;
imshow(data(:,:,1,4),[]);

hold;

disp('display d(1).imData');
figure;
imshow(d(1).imData(:,:,1,1),[]);

