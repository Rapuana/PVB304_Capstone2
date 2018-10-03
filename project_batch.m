%% Project batch file



% Run the batch 
clear, clc;
tic
job = batch('project2', 'Pool', 3, ...
    'AttachedFiles', {'rubinovmouse.mat','brainwaves_long_varyW.m'});

toc

load(job)
