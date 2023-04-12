% run recenter channels before running this script to 
% adjust the BIDS locations from the downloaded datasets

clear
bidsRepo = 'ds003061';

%% ---------------
%% End of settings
%% ---------------

if strcmpi(bidsRepo, 'ds002718')
    bidsRepo      = 'ds002718';
    bidsTask      = 'task-FaceRecognition';
    bidsSess      = '';
    bidsRun       = '';
    conditions    = { 'famous_new' 'scrambled_new' }; % types 5 and 17
    timeRanges = { [250 350] };
elseif strcmpi(bidsRepo, 'ds002680')
    bidsRepo      = 'ds002680';
    bidsTask      = 'task-gonogo';
    bidsSess      = 'ses-1';
    bidsRun       = '';
    conditions    = { 'animal_target_correct' 'animal_distractor_correct' };
    timeRanges = { [350 450] };
elseif strcmpi(bidsRepo, 'ds003061')
    bidsRepo      = 'ds003061';
    bidsTask      = 'task-P300';
    bidsSess      = '';
    bidsRun       = '2';
    conditions    = { 'oddball_with_reponse' 'standard' };
    timeRanges = { [400 500] };
else
    error('Unknown bids repo')
end
rng(now)

% add necessary path (your path might be different
if ~exist('eeg_checkset.m', 'file')
    addpath('~/eeglab');
    eeglab;
end

participants = readtable(fullfile( bidsRepo, 'participants.tsv'), 'filetype', 'delimitedtext');
nParticipants = size(participants,1);
pca_single = cell(1, nParticipants);
pca_double = cell(1, nParticipants);
mir_single = cell(1, nParticipants);
mir_double = cell(1, nParticipants);
parfor iSubject = 1:nParticipants

    subject = participants{iSubject,1}{1};
    pca_single_tmp = zeros(1,3);
    pca_double_tmp = zeros(1,3);
    mir_single_tmp = zeros(1,3);
    mir_double_tmp = zeros(1,3);
    for iRun = 1:3
        fileName = makebidsfile( '.', bidsRepo, subject, bidsSess, bidsTask, num2str(iRun));
        fprintf('Processing file %s\n', fileName);
        filePath = fileparts(fileName);
    
        EEG = pop_loadset(fileName);
        removeChans = {'EXG1','EXG2','EXG3','EXG4','EXG5','EXG6','EXG7','EXG8','GSR1','GSR2','Erg1','Erg2','Resp','Plet','Temp','061','062','063','064'};
        EEG = pop_loadset('filename',fileName);
        EEG = pop_eegfiltnew(EEG, 'locutoff',0.5);  
        EEG = pop_select( EEG, 'nochannel',removeChans); % list here channels to ignore
    
        EEG = pop_reref(EEG, []);
        [pc,eigvec,sv1] = runpca(single(EEG.data));
        [pc,eigvec,sv2] = runpca(double(EEG.data));
        pca_single_tmp(iRun) = sv1(end,end);
        pca_double_tmp(iRun) = sv2(end,end);

        [Y, W] = picard(EEG.data, 'maxiter', 500, 'mode', 'standard', 'verbose', true);
        mir_single_tmp(iRun) = getmir(double(W), double(EEG.data));

        [Y, W] = picard(double(EEG.data), 'maxiter', 500, 'mode', 'standard', 'verbose', true);
        mir_double_tmp(iRun) = getmir(double(W), double(EEG.data));
    end
    pca_single{iSubject}  = pca_single_tmp;
    pca_double{iSubject}  = pca_double_tmp;
    mir_single{iSubject}  = mir_single_tmp;
    mir_double{iSubject}  = mir_double_tmp;
end

if iscell(pca_single)
     pca_single = [ pca_single{:} ];
     pca_double = [ pca_double{:} ];
     mir_single = [ mir_single{:} ];
     mir_double = [ mir_double{:} ];
end

figure; 
hist2([pca_single], [pca_double]);
legend({'single', 'double'})

figure; 
plot([pca_single], [pca_double], 'b.');
[ypred, alpha, rsq, slope, intercept] = fastregress([pca_single], [pca_double], 1, 0);
axis equal
legend({'single', 'double'})

fprintf('Mean and std single: %1.4f (%1.4f)\n', mean([pca_single]), std([pca_single]));
fprintf('Mean and std double: %1.4f (%1.4f)\n', mean([pca_double]), std([pca_double]));
fprintf('Sign test: %g\n', signtest([pca_single]-[pca_double]));
