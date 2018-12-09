%% Project 3
clear, clc;

format compact
format short


% Go to project
directory = cd;
addpath(genpath(directory))

% Author: Samuel Dudley

% Set output directory as string
outdir = [directory filesep 'OUTPUTS'];

% Load the simulation
load rubinovmouse.mat


% Set label name 
label = 'mouse';



% Set variables for the program
% This normalises the weightings
normW = rubinovmouse.W/max(rubinovmouse.W(:));

%% Set GLOBAL parameters and perturbations
global ictime segments

ictime = 60; segments = 50;


clear t_perturb strength_perturb node 

t_perturb = 0; % 121; % integer>0 or NaN
strength_perturb = 0; % 0.05; % real or NaN


% SET NODE of interest
node = 24; % % integer 

% Set perturb here. 
perturb = false; % true;

if perturb
    t_perturb = 150; % integer>0 or NaN
    strength_perturb = 0.2; % real or NaN
end

%% Define the parameters here

delay = 1;
coupling = 0.59;


%% Now we want to run sam_perturbations_BW

% For this instance there should be no perturbation
% sam_perturbations_BW(normW, coupling, delay, label, outdir) %, ...
    %t_perturb, strength_perturb, node)

if perturb
    sam_perturbations_BW(normW, coupling, delay, label, outdir, ...
                            t_perturb, strength_perturb, node)
else
    sam_perturbations_BW(normW, coupling, delay, label, outdir)
end


%% This is where we load the data
loadedfile = loadname(coupling, delay, strength_perturb, node, 'mouse');
File = dir([outdir filesep 'Run*.*']);
ListOfFiles = {File.name};

if any(strcmp(ListOfFiles, loadedfile))
    in = load([outdir filesep loadedfile])
else
    error('No files with these variables. Try creating some')
end
    
%%

close all;
% Try node 4.
if perturb
    tp = in.time;
    yp_perturb_pos = in.soln';
else
    t = in.time;
    yp = in.soln';
end
node = in.node;

%%



figure
plot(t, yp(:,node))

hold on
% figure
plot(tp,yp_perturb(:,node))

plot(t_perturb, yp_perturb(t_perturb,node),'*r')




legend('Normal', 'Perturbed', 'time')
title(['perturbation str ' num2str(strength_perturb)])

%%


dnode = dnode + 1;
dnode = node


figure(1)
plot(t, yp(:,dnode), '-g', 'linewidth', 2)
% title([num2str(dnode), ', Non perturbed'])

rubinovmouse.name(dnode)

hold on
plot(t, yp_perturb(:, dnode), '-r', 'linewidth', 2)
% title([num2str(dnode), ', Perturbed'])


plot(t, yp_perturb_pos(:, dnode), '-b', 'linewidth', 2)

legend('Normal', 'Neg - perturbed', 'Pos - perturbed')

xlim([150 350])

title('Node 24, both positive and negative perturbations of 0.2', 'interpreter', 'latex', 'fontsize', 16)


%%
figure
plot(tp, yp_perturb(:,node))
hold on
plot(t, yp(:,node))

%%



sphereanim_plot(yp)

















%% FUNCTIONS


function name = loadname(coupling, delay, perturb, node, label)
% A function that takes inputs and creates loadname
global ictime segments

c = coupling;
if perturb
    basename = ['RunP_ictime' num2str(ictime) '_seg' num2str(segments)];
else
    basename = ['Run_ictime' num2str(ictime) '_seg' num2str(segments)];
end

if perturb
    name = ...
        sprintf('%s_d%.fms_%s_coupling%.3f_node%g_str%g.mat', ...
                basename, delay, label, c, node, perturb);
else
    name = ...
        sprintf('%s_d%.fms_%s_coupling%.2f.mat', ...
                basename, delay, label, c);
end

end