
% This function processes the GS data set
% function fitData_GS_T1 (SUBJECT_ID, T1path, indeces)

function fitData_GS_T1 (subjectID, T1path, indeces_complex, varargin)

if exist('log_T1_analysis_complex','file')
    ans0 = input(['log_T1_analysis exists. Overwrite or append [o/a]?  '],'s');
    
    if ans0=='o',
        delete('log_T1_analysis');
    end;
end;

if nargin > 3
    complex_flag = varargin{1};
else
    complex_flag = 1;
end

diary ('log_T1_analysis_complex')

%T1path = '/data/mril/mril11/stikov/Scans/T1_mapping/gu_ye_20110727_123110/GS/data_tmp';

% Where to find the DICOM data
%loadpath = [T1path 'singleslicedicomfiles/'];

%indeces = [17 19 21 23 25 27 29 31 33 35]; % Selecting only the right indeces


indeces = downsample(indeces_complex, 2);
%[T1path '/data_tmp']


loadpath = T1path;

% Where to save the .mat
%savename = [T1path 'data/' 'TestSingleSlice'];

%% Complex Fit

% Now trying to do a magnitude/phase data fit
% For that we need to use convert_mp_ri.m
%

if (complex_flag==1)
    getDataComplex(T1path, indeces_complex, subjectID);
    
    dataStr = fullfile (T1path, 'dataComplex.mat');
    savename = fullfile (T1path, 'resultsComplex.mat');
    %
    T1ScanExperiment(dataStr, savename, 'RD-NLS');
    T1FitDisplayScanLabel(dataStr, savename);
    %T1FitDisplayScan(dataStr, savename);
    %
    
else
    
    %% Magnitude Fit
    % %
    getDataMag (T1path, indeces);
    %
    dataStr = fullfile (T1path, 'data.mat');
    savename = fullfile (T1path, 'results.mat');
    
    T1ScanExperiment(dataStr, savename, 'RD-NLS-PR');
    
    % T1FitDisplayScan(dataStr, savename);
    % Call T1FitDisplayScanLabel for automated label input
   T1FitDisplayScanLabel(dataStr, savename);
    %T1FitDisplayScan(dataStr, savename);
end

diary off;

