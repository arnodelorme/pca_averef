% run recenter channels before running this script to 
% adjust the BIDS locations from the downloaded datasets

clear
bidsRepo = 'ds003061';
bidsTask      = 'task-P300';
bidsSess      = '';
icaType  = '';
parallel = 39;

try
    parpool(parallel);
catch
end

%% ---------------
%% End of settings
%% ---------------

% add necessary path (your path might be different
if ~exist('eeg_checkset.m', 'file')
    addpath('~/eeglab');
    eeglab;
end

participants = readtable(fullfile( bidsRepo, 'participants.tsv'), 'filetype', 'delimitedtext');
nParticipants = size(participants,1);
pca_single1 = zeros(1, nParticipants*3);
pca_single2 = zeros(1, nParticipants*3);
pca_double1 = zeros(1, nParticipants*3);
pca_double2 = zeros(1, nParticipants*3);
mir_single  = zeros(1, nParticipants*3);
mir_double  = zeros(1, nParticipants*3);
mir_double2  = zeros(1, nParticipants*3);
parfor iSubjectRun = 1:nParticipants*3

    iSubject = floor((iSubjectRun-1)/3)+1;
    iRun     = mod(iSubjectRun-1,3)+1;

    subject = participants{iSubject,1}{1};
    fileName = makebidsfile( '.', bidsRepo, subject, bidsSess, bidsTask, num2str(iRun));
    fprintf('Processing file %s\n', fileName);
    filePath = fileparts(fileName);

    EEG = pop_loadset(fileName);

    % preprocess data
    removeChans = {'EXG1','EXG2','EXG3','EXG4','EXG5','EXG6','EXG7','EXG8','GSR1','GSR2','Erg1','Erg2','Resp','Plet','Temp','061','062','063','064'};
    EEG = pop_loadset('filename',fileName);
    EEG = pop_eegfiltnew(EEG, 'locutoff',0.5);  
    EEG = pop_select( EEG, 'nochannel',removeChans); % list here channels to ignore

    % rereference
    EEG = pop_reref(EEG, []);
    [pc,eigvec,sv1] = runpca(single(EEG.data)); pca_single1(iSubjectRun) = sv1(end,end);
    [pc,eigvec,sv2] = runpca(double(EEG.data)); pca_double1(iSubjectRun) = sv2(end,end);
    pca_single2(iSubjectRun) = lowest_eig(single(EEG.data));
    pca_double2(iSubjectRun) = lowest_eig(double(EEG.data));

    tmpdata = runpca(double(EEG.data), EEG.nbchan-1);
    [W,S] = runica_single(single(tmpdata));
    mir_single(iSubjectRun) = getmir(double(W*S), double(tmpdata));

    [W,S] = runica(double(tmpdata));
    mir_double(iSubjectRun) = getmir(double(W*S), tmpdata);

    [W,S] = runica(double(EEG.data(1:end-1,:)));
    mir_double2(iSubjectRun) = getmir(double(W*S), double(EEG.data(1:end-1,:)));
end

printvar(pca_single1);
printvar(pca_single2);
printvar(pca_double1);
printvar(pca_double2);
printvar(mir_single);
printvar(mir_double);
printvar(mir_double2);
save('-mat', [ 'pca_' icaType '_' datestr(now, 30) '.mat'], 'pca_single1', 'pca_single2', 'pca_double1', 'pca_double2', 'mir_single', 'mir_double', 'mir_double2');


