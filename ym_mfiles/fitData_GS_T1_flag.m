
% This function processes the GS data set
% function fitData_GS_T1_flag (SUBJECT_ID, T1path, indeces_complex, complex_flag)
% Input: subjectID
%        T1path:
%        indeces_complex: the image numbers for all the T1 data (include
%        mag and phase images)
%        complex_flag: 0 for magnitude fit, 1 for complex fit
% Output: plots showing results
% Written by Stanford University
% Modified by Yuhan Ma, Jan 22, 2013
% Modifications: added the complex_flag
%

function fitData_GS_T1_flag (subjectID, T1path, indeces_complex, complex_flag)

if (~exist('complex_flag','var'))
    complex_flag = 0;
end

if exist('log_T1_analysis_complex','file')
    ans0 = input(['log_T1_analysis exists. Overwrite or append [o/a]?  '],'s');
    
    if ans0=='o',
        delete('log_T1_analysis');
    end;
end;



diary ('log_T1_analysis_complex')

%T1path = '/data/mril/mril11/stikov/Scans/T1_mapping/gu_ye_20110727_123110/GS/data_tmp';

% Where to find the DICOM data
%loadpath = [T1path 'singleslicedicomfiles/'];

%indeces = [17 19 21 23 25 27 29 31 33 35]; % Selecting only the right indeces



%[T1path '/data_tmp']


loadpath = T1path;

% Where to save the .mat
%savename = [T1path 'data/' 'TestSingleSlice'];
%% Magnitude Fit
if complex_flag == 0
    indeces = downsample(indeces_complex, 2);
    
    getDataMag (T1path, indeces);
    
    dataStr = fullfile (T1path, 'data.mat');
    savename = fullfile (T1path, 'results.mat');
    
    T1ScanExperiment(dataStr, savename, 'RD-NLS-PR');
    
    % T1FitDisplayScan(dataStr, savename);
    % Call T1FitDisplayScanLabel for automated label input
    T1FitDisplayScanLabel(dataStr, savename);
    
    
    
    %% Complex Fit
elseif complex_flag == 1
    
    % Now trying to do a magnitude/phase data fit
    % For that we need to use convert_mp_ri.m
    
    getDataComplex(T1path, indeces_complex, subjectID);
    
    dataStr = fullfile (T1path, 'dataComplex.mat');
    savename = fullfile (T1path, 'resultsComplex.mat');
    %
    T1ScanExperiment(dataStr, savename, 'RD-NLS');
    T1FitDisplayScanLabel(dataStr, savename);
    
else
    fprintf('Invalid complex_flag values. Please use 0 for magnitude fit, 1 for complex fit\n');
end


diary off;

