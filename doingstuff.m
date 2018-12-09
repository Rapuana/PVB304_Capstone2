%% OBSERVING THE AFFECT OF PERTURBAITON ON OTHER NODES

global ictime segments delays label c strength_perturb outdir t y
% Load the data
load('baseline_59coupling_1delay.mat')

t = baseline.time';
y = baseline.soln';

normW = baseline.normW;
c = baseline.coupling;
delays = baseline.delays;

label = baseline.label;

strength_perturb = 0.2;

ictime = baseline.ictime;
segments = baseline.segments;


clear baseline

%% This is the section to produce the figures

figure(1)
plot(t,y)
xlabel('\textbf{Time} (ms)', 'Interpreter', 'latex', 'fontsize', 14)
ylabel('\textbf{Voltage} (mV)', 'interpreter', 'latex', 'fontsize', 14)

axis([1500 2000 -0.45 0.2])

hold on
% Can we plot the average?

yavg = mean(y,2);

% figure(2)
plot(t,yavg, '-r', 'linewidth', 2)
figure(1)


% Can we figure out which ones fire above all the rest?


threshold = 0;
plot([1500 2000], [threshold threshold], '-k', 'linewidth', 2)

%%

% Let's plot it one node at a time

nd = 11; % Node of interest

trange = 500:length(t);


figure(2)
plot(t(trange), y(trange, nd))


%% ANALYSIS

[~, locs] = findpeaks(y(trange,nd));

cycles = diff(locs);
meanCycle = mean(cycles);

meanCycle = 1:112;

for ii = 1:112
    [~, locs] = findpeaks(y(trange,ii));
    meanCycle(ii) = mean(diff(locs));
end
    


% Find spikes (we define spikes to be greater than the threshold we set.
%% LOOP it Br\"other 
spikehist_baseline = [];

figure(3)
hold on

for ii = 1:112 % All the nodes we have
    [pks,locs] = findpeaks(y(trange,ii), 'minpeakheight', threshold);
    cycles = diff(locs);
    if ~isempty(pks)
        jj = find(max(pks));
        spikehist_baseline = [spikehist_baseline; [ii pks(jj) locs(jj)]];
        plot(t(trange),y(trange,ii))
    end
end




%% let's have a look
figure(4)
plot(t(trange), y(trange,14))

ii = find(spikehist_baseline(:,1) == 14);

% This places the perttime at roughly 5 cycles before the peaks spike.
% What do we want to hopefully observe. Does this spike occur again?
% If so when?

pert_time = spikehist_baseline(ii,3) - floor(5*meanCycle(14)) + trange(1);




%% Make the big boi loop adding perturbations at the required spots

directory = cd;
outdir = [directory filesep 'OUTPUTS' filesep 'OUTPUTS_PERTURBED'];


if ~true

for jj = 1:length(spikehist_baseline)
    pert_time = spikehist_baseline(jj,3)-floor(5*meanCycle(14))+trange(1);
    sam_perturbations_BW(normW,c,delays,label, outdir,...
                      pert_time,strength_perturb,spikehist_baseline(jj,1))
end
end

%% Now we want to compare these with the og

% Fixed params

% We only change node
node = 6;


% Load values
basename = ['RunP_ictime' num2str(ictime) '_seg' num2str(segments)];
name = ...
    sprintf('%s_d%.fms_%s_coupling%.3f_node%g_str%g.mat', ...
    basename, delays, label, c, node, strength_perturb);

input = load([outdir filesep name]);

y6 = input.soln';


t_perturb = input.t_perturb;




%%
figure
plot(t(t_perturb:end), y(t_perturb:end,node))
hold on
plot(t(t_perturb:end), y6(t_perturb:end,node))

legend('normal', 'perturbed')
%%
figure(1)
plot(t(t_perturb:end), y(t_perturb:end,:))
hold on
plot(t(t_perturb:end), y(t_perturb:end,node), '-k', 'linewidth', 1.5)

figure(2)
plot(t(t_perturb:end), y6(t_perturb:end,:))
hold on
plot(t(t_perturb:end), y6(t_perturb:end,node), '-k', 'linewidth', 1.5)

%% What is that weird one?


%% What about the rest? ~ Worry later lol
figure
plot(t(trange), y(trange, 6))

%% Let's create something to automatically figure out the stuff. namely what is the lowest
numArrays = 112;
yp = cell(numArrays,1);
tp = cell(numArrays,1);
for j = 1:length(spikehist_baseline)
    jj = spikehist_baseline(j,1);
    
    input = loading(jj);
    yp{jj} = input.soln';
    tp{jj} = input.t_perturb;
end

%%

% Now we need to look over the entire function and figure out which ones
% give us the worst values

min_nodes = 1:length(spikehist_baseline);

for i = 1:length(min_nodes)
    ii = spikehist_baseline(i,1);
    [~,I] = min(min(yp{ii}(500:end,:)));
    min_nodes(i) = I;
end




%%
figure
plot(t, yp{7})

figure
plot(t(600:800), yp{7}(600:800, 7))
hold on
plot(t(600:800), y(600:800, 7))
legend('perturbed', 'normal')
tp{7}


%% 
for ii = 1:length(spikehist_baseline)
    i = spikehist_baseline(ii,1);
    figure
    plot(t(tp{i}-50:tp{i}+100), y(tp{i}-50:tp{i}+100,i))
    hold on
    plot(t(tp{i}-50:tp{i}+100), yp{i}(tp{i}-50:tp{i}+100,i))
    plot(t(tp{i}), y(tp{i},i), '*r')
    legend('normal', 'perturbed', 'perturbation')
    title(num2str(i))
end

%%
nd = 94;
global nd
figure
plot(t(tp{nd}-50:tp{nd}+100), y(tp{nd}-50:tp{nd}+100,nd))
    hold on
    plot(t(tp{nd}-50:tp{nd}+100), yp{nd}(tp{nd}-50:tp{nd}+100,nd))
    plot(t(tp{nd}), y(tp{nd},nd), '*r')
    legend('normal', 'perturbed', 'perturbation')
    title(num2str(nd))
    
    
%% Now we want to add perturbations of different kinds. Let's add a perturbation at 1670

t_perturb = 1670;
strength_perturb = 0.05;
if ~true
sam_perturbations_BW(normW,c,delays,label, outdir,...
                      t_perturb,strength_perturb,nd)
end
%%
input = loading(nd);


ytest = input.soln';

%% now plot 
figure(3)
plot(t(tp{nd}-50:tp{nd}+100), y(tp{nd}-50:tp{nd}+100,nd))
    hold on
    plot(t(tp{nd}-50:tp{nd}+100), ytest(tp{nd}-50:tp{nd}+100,nd))
    plot(t(t_perturb), y(t_perturb,nd), '*r')
    legend('normal', 'perturbed', 'perturbation')
    title(num2str(nd))
    
figure(1)
plot(t, ytest)
hold on 
plot(t, ytest(:,nd), '-k', 'linewidth', 1.5)
xlim([1600 2000])
title('strength 0.05')


figure(2)

plot(t, y)
hold on
plot(t, y(:,nd), '-k', 'linewidth', 1.5)
xlim([1600 2000])
ylim([-0.5 0.4])


%% Try again with higher strength
t_perturb = 1670;
strength_perturb = 0.05;

sam_perturbations_BW(normW,c,delays,label, outdir,...
                      t_perturb,strength_perturb,nd)
                  
                  
                  
%%
input = loading(nd);


ytest20 = input.soln';


figure(3)
plot(t(tp{nd}-50:tp{nd}+100), y(tp{nd}-50:tp{nd}+100,nd))
    hold on
    plot(t(tp{nd}-50:tp{nd}+100), ytest20(tp{nd}-50:tp{nd}+100,nd))
    plot(t(tp{nd}-50:tp{nd}+100), ytest(tp{nd}-50:tp{nd}+100,nd))
    plot(t(t_perturb), y(t_perturb,nd), '*r')
    legend('normal', 'perturbed 0.20', 'perturbed 0.05', 'perturbation')
    title(num2str(nd))
    
figure(4)
plot(t, ytest20)
hold on 
plot(t, ytest20(:,nd), '-k', 'linewidth', 1.5)
plot(t, ytest(:,nd), '-r', 'linewidth', 1.5)
xlim([1600 2000])
title('strength of 0.20')

figure(2)
plot(t, y)
hold on
plot(t, y(:,nd), '-k', 'linewidth', 1.5)
xlim([1600 2000])
ylim([-0.5 0.4])



%% Interesting phenomena
find(ytest(2920,:) > -0.05)

xpks = [];
xlocs = [];
xnode = [];
for i = 1:112
    [pks, locs] = findpeaks(y(2900:end,i), 'minpeakheight', -0.05)
    xpks = [xpks; pks];
    xlocs = [xlocs; locs];
    mynode = i .* ones(length(locs),1);
    xnode = [xnode; mynode];
end

xpks20 = [];
xlocs20 = [];
xnode20 = [];
for i = 1:112
    [pks, locs] = findpeaks(ytest20(2900:end,i), 'minpeakheight', -0.05)
    xpks20 = [xpks20; pks];
    xlocs20 = [xlocs20; locs];
    mynode = i .* ones(length(locs),1);
    xnode20 = [xnode20; mynode];
end

xpks05 = [];
xlocs05 = [];
xnode05 = [];
for i = 1:112
    [pks, locs] = findpeaks(ytest(2900:end,i), 'minpeakheight', -0.05)
    xpks05 = [xpks05; pks];
    xlocs05 = [xlocs05; locs];
    mynode = i .* ones(length(locs),1);
    xnode05 = [xnode05; mynode];
end


%%
figure(7)
plot(t(2850:end), y(2850:end, [4 9 11 14 65]));
title('normal')
figure(8)
plot(t(2850:end), ytest(2850:end, [4 9 11 14 65]));
title('perturbation of 0.05')
figure(9)
plot(t(2850:end), ytest20(2850:end, [4 9 11 14 65]));
title('perturbation of 0.20')

%%

figure(6)
hold on
for i = 1:length(xnode)
    figure(6)
    
    plot(t(2850:end), y(2850:end, xnode(i)))
end



%%
nod = find((ytest20(1806,:)) > -0.04)

figure
plot(t, ytest20(:,nod))
hold on
plot(t, y(:,nod))
xlim([1600 3000])



%% a different strength and different location?

nd = 94;

t_perturb = 1690;
strength_perturb = -0.05;

sam_perturbations_BW(normW,c,delays,label, outdir,...
                      t_perturb,strength_perturb,nd)
                  
%%
input = loading(nd);


ytestn05 = input.soln';


figure(3)
plot(t(tp{nd}-50:tp{nd}+100), y(tp{nd}-50:tp{nd}+100,nod))
    hold on
    plot(t(tp{nd}-50:tp{nd}+100), ytest(tp{nd}-50:tp{nd}+100,nod))
    plot(t(tp{nd}-50:tp{nd}+100), ytestn05(tp{nd}-50:tp{nd}+100,nod))
    plot(t(t_perturb), y(t_perturb,nod), '*r')
    legend('normal', 'perturbed 0.20', 'perturbed -0.05', 'perturbation')
    title(num2str(nd))

    
    
figure(5)
plot(t(tp{nd}-50:tp{nd}+100), ytestn05(tp{nd}-50:tp{nd}+100,:))
figure(6)
plot(t(tp{nd}-50:tp{nd}+100), y(tp{nd}-50:tp{nd}+100,:))

%% 
nod = find(ytestn05(1698,:) > 0)
%% functions
directory = cd;
outdir = [directory filesep 'OUTPUTS' filesep 'OUTPUTS_PERTURBED'];

% A function
function input = loading(node)
global ictime segments delays label c strength_perturb outdir 
% Load values
basename = ['RunP_ictime' num2str(ictime) '_seg' num2str(segments)];
name = ...
    sprintf('%s_d%.fms_%s_coupling%.3f_node%g_str%g.mat', ...
    basename, delays, label, c, node, strength_perturb);
input = load([outdir filesep name]);
end

function figures(yx, tperturb,node)
global t y 

figure(1)
plot(t(tperturb:end), y(tperturb:end,:))
hold on
plot(t(tperturb:end), y(tperturb:end,node), '-k', 'linewidth', 1.5)
title('No Perturbation')

figure(2)
plot(t(tperturb:end), yx(tperturb:end,:))
hold on
plot(t(tperturb:end), yx(tperturb:end,node), '-k', 'linewidth', 1.5)
title(['Perturbation at ' num2str(tperturb) ' ms for node ' num2str(node)])
end




%


