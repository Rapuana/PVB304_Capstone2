%% Project v2
clear, clc;

format compact
format short

% Go to project
cd N:\Project
addpath(genpath('N:\Project'))

% Author: Samuel Dudley

% Set output directory as string
directory = 'C:\Users\samD\Documents\project\Final_Outputs_Project1';

% Load the simulation
load rubinovmouse.mat


% Set variables for the program
normW = rubinovmouse.W/max(rubinovmouse.W(:));

%% Define the parameters here

delays = 1;
couplings = 0.61:0.01:0.69;

% overwrite if necessary



%% Create the data here
% For some reason this really hates a coupling of 1 with some delay that is
% close to the coupling value. 

% strangely enough it is an index error (positive of logical values).


% HISTORY of failures

% Failed on 
% delay = 1.0, coupling = 1.0  --  x9
% delay = 0.9, coupling = 1.0  --  x8

% redo values
d = 1; c = 0.7;

% brainwaves_long_varyW(normW, c, d, 'mouse', directory);

for i = 1:length(delays)
    for k = 1:length(couplings)
        brainwaves_long_varyW(normW, couplings(k), delays(i), 'mouse', directory);
    end
end



% LOAD VALUES
%         c,   d;
POI =  [0.6, 1.0;
        0.7, 0.7;
        0.7, 1.0;
        0.7, 0.8;
        1.0, 2.0;
        0.6, 0.4];

% for i = 1:length(POI)
    
% What values do we want?
%%

close all
couplings = [0.58 0.60 0.54 0.64];
i = 1:length(couplings)

    coupling = 0.64;
    delay = 1;
    
runimages = false; % Do we want to produce and save images?
runvideos = false; % Do we want to produce and save videos?
trange = [  1300
            1400]; % For images/videos
    
% coupling = 0.65; delay = 1;
for i = 1 % HIDE ALL THIS FILE READING STUFF
% name variables to create name
ictime = 50; % [500, 50] long/short
segments = 50;
scanlag = delay;
label = 'mouse';
c = coupling;
trials = 1;
outdt = 1;

% All code below is to load the file with the parameters above
% write file name
basename=['Wtesting_ictime' num2str(ictime) '_seg' num2str(segments) '_outdt' num2str(outdt)];
loadedfile = sprintf('%s_d%.8fms_%s_coupling%.4f_trial%i.mat', basename, scanlag, label, c,trials);

% load list of all files -> check one exists if not, make one?
File = dir([directory, '\Wtesting*.*']);
ListOfFiles = {File.name};
prompt = 'Output of variables does not exist, would you like to create one? Y/N';

if any(strcmp(ListOfFiles, loadedfile))
    in = load([directory, '\', loadedfile])
else
    runcode = input(prompt,'s');
    if runcode == 'Y' || runcode == 'y'
        brainwaves_long_varyW(normW, coupling, delay, 'mouse', directory);
        in = load([directory, '\', loadedfile])
    elseif runcode == 'N' || runcode == 'n'
        sprintf('Coupling = %.2f, Delay = %.2f.', coupling, delay)
        error('File does not exist for chosen coupling and delay values')
    else
        error('Please enter valid input of Y/N or try changing Coupling and Delay')
    end
end

% Show the data
% Name Scheme

long = NaN;

if ictime == 50
    long = '_ic50';
elseif ictime == 500
    long = '_long';
else
    long = ['_ict_', num2str(ictime)];
end

fps = 10;
IMGname = sprintf('TimeSeries_coupling_%.2f_delay_%.2f', coupling, delay);
vidname = sprintf('video_coupling_%.2f_delay_%.2f_fps_%.0f', coupling, delay, fps);
end
clc
close all
t64 = in.time;
yp64 = in.soln';
%


if runimages
figure
plot(t,yp)
xlabel('\textbf{Time} (ms)', 'Interpreter', 'latex')
ylabel('\textbf{Voltage} (mV)', 'interpreter', 'latex')

yprange = [min(min(yp)), max(max(yp))]; 
axis([min(t) 1.1*max(t) yprange(1) yprange(2)])
y_pos = yprange(2) - (yprange(2)-yprange(1)) * 0.08;
x_pos = 0.7*max(t);
label_str = {['\textbf{Coupling} = ', num2str(coupling)]; ['\textbf{Delay} = ', num2str(delay), ' ms']};
text(x_pos, y_pos, label_str, 'Interpreter', 'latex', 'fontsize', 12)
%
saveas(gcf, [directory, '\', IMGname, long, '.png'])

% Excerpt Image


time = [1300, 2100];
trange = time;
% close all

figure
plot(t,yp)
xlabel('\textbf{Time} (ms)', 'Interpreter', 'latex')
ylabel('\textbf{Voltage} (mV)', 'interpreter', 'latex')

yprange = [min(min(yp(trange(1):trange(2),:))), max(max(yp(trange(1):trange(2),:)))]; 
axis([trange(1) trange(1)+1.1*(trange(2)-trange(1)) yprange(1) yprange(2)])
y_pos = yprange(2) - (yprange(2)-yprange(1)) * 0.08;
x_pos = 0.7*diff(time)+time(1);
label_str = {['\textbf{Coupling} = ', num2str(coupling)]; ['\textbf{Delay} = ', num2str(delay), ' ms']};
text(x_pos, y_pos, label_str, 'Interpreter', 'latex', 'fontsize', 12)

saveas(gcf, ['C:\Users\samD\desktop', '\', IMGname, '_ROI.png'])
end
% Videos
if runvideos

% coupling, delay
% sphereanim_plot(yp, 0.02, 'spring')


%
az = 0; el = 0; roll = -90;

% Save the video
% brainwaves_vid_spheres([directory, '\', vidname, long, '.mp4'],yp, fps, 'jet', az, el, roll)

% end



% %% image excerpt


% ROI refers to 1300 to 2100 ms and is the focus of what I am doing now

time = [1300 2100];

% %%
fps = 10;
vidname1 = sprintf('EXTRACT_video_coupling_%.2f_delay_%.2f_fps_%.0f_t1_%.f_length_%.f', coupling, delay, fps, time(1), diff(time));

%brainwaves_vid_spheres([directory, '\', vidname1, '_sideview.mp4'],yp, fps, 'jet', 0, -90, 0, time)
brainwaves_vid_spheres(['C:\Users\samD\desktop', '\', vidname1, '.mp4'], yp, fps, 'jet', az, el, roll, time);
% for i = 1:length(t)
%     sphere_picture(yp, time(i),'jet')
% end


%
% sphere_picture(yp,1460, 'jet')
end
%% 
% What do we want to do?

% We want to sort the mouse connectomes such that it is ordered from 'y'.

% thus using 
[~, I] = sort(rubinovmouse.XYZ(:, 1));
% We have the values in sorted order, Y
% and the index as to where they were.

t1 = 1830;
t2 = 1900;

%%


d = yp64(t1:t2,:);

a = zeros(size(d));

for i = 1:length(yp(1,:))
    a(:,i) = d(:,I(i));
end

a64 = a;
%%
subplot(5, 1, 5)
colormap jet
imagesc([t1 t2], [1 112], a64', [-0.3 0]) % Ordered by 'X' positioning (useful for determining waves)
caxis = [-0.3 0]

titlename = sprintf('c = %.2f', 0.64);
title(titlename, 'interpreter', 'latex')


%%
cbar = colorbar;
ylabel(cbar, 'V')

ylabel('\textbf{Connectome}', 'interpreter', 'Latex', 'FontSize', 16)
xlabel('\textbf{Time} (ms)', 'interpreter', 'latex', 'FontSize', 16)



% figure
% imagesc(d') % Original, untouched 
% colormap jet

% 
% figure 
% surf(a')
% 
% ylabel('\textbf{Connectome}', 'interpreter', 'Latex', 'FontSize', 16)
% xlabel('\textbf{Time} (ms)', 'interpreter', 'latex', 'FontSize', 16)
% zlabel('\textbf{$\mu_e$} (V)', 'interpreter', 'latex', 'Fontsize', 16)

%%
close all
figure
imagesc(log10(normW))
cbar = colorbar
ylabel(cbar, '$\log(c)$', 'interpreter', 'latex', 'fontsize', 16)

saveas(gcf, 'C:/users/samD/desktop/Connectivity_Matrix.png')

%% Save the images here

% figure, plot(sum(normW))
% figure, imagesc(log10(normW))
% figure, plot(sum(normW>0))
% figure, figure, plot(yp), hold on, plot(yp(:,[9 12]),'linewidth',5)
% figure, plot(yp(end,:))
% rubinovmouse
% 
% rubinovmouse.name{1}
% ans =
%     'Ammon's horn'
% rubinovmouse.name{2}
% ans =
%     'Dentate gyrus'
% figure, imagesc(corr(yp(1500:end,:)))
% ysnippet=yp(1500:end,:);
% yhilb=hilbert(ysnippet-mean(yysnippet));
% Undefined function or variable 'yysnippet'. 
% Did you mean:
% yhilb=hilbert(ysnippet-mean(ysnippet));
% figure, plot(ysnippet(:,1)), hold on, plot(abs(yhilb(:,1)))
% figure, plot(ysnippet(:,1)), hold on, plot(abs(yhilb(:,1))-mean(ysnippet(:,1)))
% figure, plot(ysnippet(:,1)), hold on, plot(abs(yhilb(:,1))+mean(ysnippet(:,1)))
% figure, plot(ysnippet(:,1)), hold on, plot(angle(yhilb(:,1)))
% yphase=angle(yhilb);
% abs(mean(exp(1i*yphase)))
% ans =
%   Columns 1 through 21
%     0.0447    0.0525    0.0051    0.0072    0.0126    0.0453    0.0439    0.0160    0.0017    0.0067    0.0022    0.0055    0.0040    0.0067    0.0079    0.0199    0.0018    0.0079    0.0091    0.0010    0.0039
%   Columns 22 through 42
%     0.0052    0.0026    0.0026    0.0080    0.0208    0.0010    0.0051    0.0093    0.0064    0.0022    0.0509    0.0232    0.0537    0.0219    0.1130    0.0013    0.0033    0.0028    0.0154    0.0044    0.0006
%   Columns 43 through 63
%     0.0882    0.0053    0.0005    0.0035    0.0008    0.0060    0.0054    0.0169    0.0230    0.0129    0.0044    0.0148    0.0063    0.0008    0.0450    0.0523    0.0017    0.0104    0.0160    0.0430    0.0410
%   Columns 64 through 84
%     0.0194    0.0048    0.0035    0.0054    0.0041    0.0065    0.0038    0.0108    0.0177    0.0044    0.0103    0.0069    0.0039    0.0014    0.0032    0.0049    0.0007    0.0056    0.0185    0.0026    0.0073
%   Columns 85 through 105
%     0.0065    0.0097    0.0046    0.0480    0.0248    0.0515    0.0196    0.1134    0.0022    0.0055    0.0017    0.0127    0.0070    0.0019    0.0886    0.0027    0.0020    0.0025    0.0021    0.0037    0.0031
%   Columns 106 through 112
%     0.0142    0.0256    0.0162    0.0073    0.0178    0.0100    0.0034
% figure, plot(abs(mean(exp(1i*yphase))))
% figure, plot(abs(mean(exp(1i*yphase))))
% figure, plot(abs(mean(exp(1i*yphase),2)))
% ylim([0 1])




