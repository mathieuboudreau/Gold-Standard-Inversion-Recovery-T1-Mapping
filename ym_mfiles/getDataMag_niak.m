% getData.m
%
% Generates the .mat from the DICOM images for GS T1 mapping. 
%
% written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%  (c) Board of Trustees, Leland Stanford Junior University    
%
% Modifications: 
%               July 2013 - Mathieu Boudreau: Adapted code to use Niak
%               input.
%               

function getDataMag_niak(subjectID,irID)

% Where to save the .mat 
savename = fullfile ('data.mat');

for ii=1:length(irID)
    [d_hdr{ii},d{ii}] = niak_read_minc([subjectID '_' num2str(irID(ii)) '_mri.mnc']);
end


nbrow = size(d{1},1);
nbcol = size(d{1},2);
nbslice = size(d{1}, 3);
nbseries = length(d);

data = zeros(nbcol,nbrow,nbslice,nbseries); % Complex data
extra.tVec = zeros(1,nbseries); % One series corresponds to one TI


for k = 1:nbseries
	dataTmp = d{k};
	dataTmp = double(squeeze(dataTmp));	
	for ss = 1:nbslice
        % ym:
        % sizedata = size(data)
		data(:,:,ss,k) = dataTmp(:,:,ss)'; 
    end
	extra.tVec(k) = 1000*cell2mat(d_hdr{k}.details.acquisition.attvalue(19));
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
