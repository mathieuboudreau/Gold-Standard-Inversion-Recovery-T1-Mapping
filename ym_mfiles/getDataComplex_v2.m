% getDataComplex.m
%
% Generates the complex .mat from the magnitude/phase MINC images.
% Example:
% getDataComplex('/data/mril/mril11/stikov/Scans/stikov_nikola_20110414_161
% 304/GS_2_complex/', [12 14 16 18], 'stikov_nikola_20110414_161304');
%
% written by N. Stikov, April 2011
% Modifications by Yuhan Ma, 2012
% Modified the convert_mp_ri (saved as convert_mp_ri_ym) and getallimages (getallimages_ym) to take data with no
% time frames
% Modifications by Yuhan Ma, April 2013
% Modifications: phase data are considered

function getDataComplex_v2(T1path, indeces, subjectID);

% clear all
% close all
%
% T1path = '/data/mril/mril11/stikov/Scans/christine_ir_t1';

% Where to find the DICOM data
%loadpath = [T1path 'singleslicedicomfiles/'];

loadpath = T1path;

% Where to save the .mat
%savename = [T1path 'data/' 'TestSingleSlice'];

s = mincLoadAllSeries(T1path);
s = s(downsample(indeces, 2));

for ii=1:numel(downsample(indeces, 2))
    
    extra.tVec(ii) = s(ii).inversionTime;
    
end

savename = fullfile (T1path, 'dataComplex.mat');

currpath = pwd;
cd(loadpath)
%addpath([T1path '/mfiles']);


% First we need to create a numframes dimension, so that convert_mp_ri
% works
%
% for ii=indeces
%     inputFile = [subjectID '_' num2str(ii) '_mri.mnc.gz'];
%     outputFile = [subjectID '_' num2str(ii) '_mri.mnc'];
%
%     eval (['!mincconcat -concat_dimension time ' inputFile ' ' outputFile]);
% end

% ym: I modified the convert_mp_ri (saved as convert_mp_ri_ym) and getallimages (getallimages_ym) to take data with no
% time frames so the above codes are commented out.

% Then we convert all the magnitude phase images into real and imaginary

jj = 1;
for ii=downsample(indeces, 2)
    magFile = [subjectID '_' num2str(ii) '_mri.mnc.gz'];
    phaseFile = [subjectID '_' num2str(ii+1) '_mri.mnc.gz'];
    
    %phase_data(:,:,1,jj) = run_prelude(subjectID, ii);
 
    [flag, realFile, imagFile] = convert_mp_ri_ym(magFile, phaseFile, 8188);
    
    realDataHandle = openimage(realFile);
    realData = getimages(realDataHandle, 1, 1);
    
    % ym: get information about height and width
    if (ii==indeces(1))
        height = getimageinfo(realDataHandle,'ImageHeight');
        width = getimageinfo(realDataHandle,'ImageWidth');
        size(realData)
    end
    closeimage(realDataHandle);
    
    imagDataHandle = openimage(imagFile);
    imagData = getimages(imagDataHandle, 1, 1);
    closeimage(imagDataHandle);
    
    %dataComplex(:,:,1,(ii - indeces(1))/2+1) = realData + 1i*imagData;
    
    % Why do we fix the slice dimension here?
    dataComplex(:,:,1,jj) = realData + 1i*imagData;
    
    jj=jj+1;
    
    
    
end


% data = reshape(dataComplex, 88, 128, 1, numel(downsample(indeces, 2)));

data = reshape(dataComplex, width, height, 1, numel(downsample(indeces, 2)));

%phase_data = angle(data);

%extra.tVec = [50 400 1100 2300 5970];


extra.T1Vec = 1:5000; % this range can be reduced if a priori information is available

% This is where you would correct the phase of certain datapoints.
% Indeed, the chopping procedure during Prescan can give an additional 180?
% phase on the whole image for certain TIs

% edited by Yuhan Ma

% num_TI = length(extra.tVec);
% 
% clims = [-3.5, 3.5];
% figure,
% for ii = 2:num_TI
%     phase_diff = phase_data(:,:,1,ii) - phase_data(:,:,1,ii-1);
%     subplot(2,3,ii-1), imagesc(squeeze(phase_diff),clims);
%     title(['Phase difference b/w TI #', num2str(ii), ' and #', num2str(ii-1)]); colorbar;hold on;
%     polarity_change = abs(phase_diff) > 3*pi/4;
% %     for jj = 1:(ii-1)
% %         data_temp = data(:,:,1,jj);
% %         data_temp(polarity_change) = -data_temp(polarity_change);
% %         data(:,:,1,jj) = data_temp;
% %     end
% end

%data(:,:,1,1:4)=-data(:,:,1,1:4);

%data(:,:,1,3:5) = -data(:,:,1,3:5);

TI = extra.tVec
save('dataComplex.mat','data','extra')
