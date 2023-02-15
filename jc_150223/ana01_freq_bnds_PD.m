    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               Get tremor boundries with Mike X Cohens approach
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data: UC San Diego Resting State EEG Data from Patients with Parkinson's Disease, OpenNeuro, 2020 (v. 1.0.3)
% Author: Julius Welzel, j.welzel@neurologie.uni-kiel.de
% EEGLab is required for this (Here EEGLab 2019_1)

clc; clear all; close all
%% Set envir

PATHJC      = fullfile('C:\Users\User\Desktop\kiel\dbs_fus_jc\202110_cohen_jw\')

addpath(genpath(PATHJC))
addpath(genpath('C:\Users\User\Documents\MATLAB\toolboxes\eeglab2022_1'))
eeglab
%% EEG analysis for 2 PD patients

%creating variables for ICA
cfg.ICA.LP = 35;
cfg.ICA.HP = 0.1;
cfg.ICA.PRUNE = 3;
cfg.ICA.resam = 250;

nms_pd = {'pd5','pd13'};

fois{1} = [13 35];

for i = 1:numel(nms_pd)
    tmp_path = fullfile(PATHJC, ['gedBounds_mxk\exsample_dataste\sub-' nms_pd{i} '\ses-off\eeg\sub-' nms_pd{i} '_ses-off_task-rest_eeg.bdf']);
    EEG = pop_biosig(tmp_path,'ref',32);

    %low-pass filter
    EEG = pop_eegfiltnew(EEG, [],cfg.ICA.LP,[],0,[],0);
    %resample to 250 Hz
    EEG = pop_resample(EEG, cfg.ICA.resam);
    %high-pass filter
    EEG = pop_eegfiltnew(EEG, [],cfg.ICA.HP,[],1,[],0);

    f = figure;pop_spectopo(EEG, 1, [0  length(EEG.data)], 'EEG' , 'percent', 100, 'freqrange',[1 40],'electrodes','off');
    
    freqbands = getBounds_JW(EEG.data,EEG.srate,[5 35],fois)

end