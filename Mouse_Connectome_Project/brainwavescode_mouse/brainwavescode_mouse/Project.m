%% First play around on a computer that works

% Samuel Dudley
dir('~/Mouse_Connectome_Project/brainwavescode_mouse/brainwavescode_mouse')

% run a simulation
load rubinovmouse.mat


% brainwaves_long_varyW(normW,coupling,delays,label,outdir)

normW = rubinovmouse.W/max(rubinovmouse.W(:));
coupling = 0.7;
delays = 1.5;


strname = brainwaves_long_varyW(normW,coupling,delays,'mouse','.')
% it will save a bunch of temporary output files, then clean them up at the end

% load in the output
in=load('Wtesting_ictime50_seg799_outdt1_d1ms_mouse_coupling0.6_trial1.mat')
t=in.time;
yp=in.soln';
%%

close all
% plot all the time series
figure, plot(t,yp)

% animate it (it plays a bit fast, might want to put a pause in there to slow it down)
sphereanim_plot(yp, 0.02, 'spring')

%%
clc
outfilename='test.mp4';
 
% vid content
frameinds=10000:10000+100-1;
yvid=in.soln(:,frameinds);
 
% setup output video
fps=10; % video frame rate
 
brainwaves_vid_spheres(outfilename,yvid,fps)