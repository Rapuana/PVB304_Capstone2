function mysavename = Perturb_brainwaves_long_varyW(normW,coupling,delays,label,outdir, perturb, t_perturb, t_length, strength_perturb, node)
% Perturb is a logical value
% t_perturb is the time in which the perturbation will occur. This is
% limited to the time steps, and must be < ictime*segments. units of ms.
% t_length is also of units ms and is the length of which the perturbation
% will occur for. We will assume that as default t_length = 1.0
% strength_perturb refers to the strength of the perturbation (mV ??? )

% global variables (shared with dde solver)
global V1 V2 V3 V4 V5 V6 V7 gCa gK gL VK VL VCa I b ani aei aie aee phi V8 V9 gNa VNa ane nse rnmda N CM vs c k_in numddevars myrand

numddevars = 3;  %number of state variables describing each unit in the network

%load connectivity matrix
%load('normalizedConnectivities.mat');

%normW=threshold_proportional(normW,thr); % could threshold using BCT if desired
%normW=double(normW>0);

%delays = [2:16];   % *********  DELAYS   ***********
%if nargin<1
%    coupling=0.6;  % ********* COUPLINGS ***********
%end
% NOTE: continues solutions in c direction

ictime=50;     %initial condition duration
segments=50; % integration time after IC transient will be segments*ictime.
             % integration is split into segments because dde23 is
             % inefficient with memory

outdt = 1;      % output timestep (in ms)


basename=['Wtesting_ictime' num2str(ictime) '_seg' num2str(segments) '_outdt' num2str(outdt)]; % output filename will start with this
%basename='test';
%label = 'W'; % goes into the output filename
if nargin<5
    outdir='~/data/networks/brainwaves/Wtesting/'; % output directory
end
if ~exist(outdir,'dir'), mkdir(outdir); end
connectivity=normW;


%ddesolver=@(varargin) dde23_JR(varargin{:}); % wrapper, in case one wishes to change dde solver
if delays==0
    ddesolver=@(varargin) dde23(varargin{:}); % wrapper, in case one wishes to change dde solver
else
    ddesolver=@(varargin) dde23(varargin{:}); % wrapper, in case one wishes to change dde solver
end

%preload brain node locations
%load('loc_AP.mat')

ntrials=1; % number of trials. 
% no longer true: trials>=2 will have perturbations as per below





% What is this ??





%For perturbations. For no perturbations set option_node=0 & rand_length=0;
option_node=0;          %=1 for random node; =0 for random time from ICs
perturbed=0;         %magnitude of perturbation
% options for option_node=0
offset=0;             %If you want to offset the time-to-perturbation
rand_length=0;        %max random initial length - perturb at offset+rand*rand_length
non_random_nde=0;     %node index
% options for option_node=1
nds=1;                  %how many nodes to perturb

rand_ics=1;    %Random initial conditions? (1=yes)
%load icssol    %Or presaved ("ics_sol")



CM = sparse(connectivity);
N = size(connectivity,1);
% set out-strength (= out-degree for binary matrices)
k_in = sum(CM)';


% MODEL PARAMS =====================================
% set model parameters
V1 = -0.01; V2 = 0.15; V3 = 0; V4 = 0.3; V5 = 0; V7 = 0; V9 = 0.3; V8 = 0.15;
gCa = 1; gK = 2.0; gL = 0.5; gNa = 6.7;
VK = -0.7; VL = -0.5; I = 0.3; b = 0.1; phi = 0.7; VNa = 0.53; VCa = 1;
ani = 0.4; vs = 1; aei = 2; aie = 2; aee = 0.36; ane = 1; rnmda = 0.25;


% more parameters: noise, coupling, modulation
nse = 0;


modn = 0;       % ********* Random parameter variance ***********
if (modn==0)
    V6 = 0.65;
else
    V6 = ones(N,1).*0.65 + modn*(rand(N,1)-0.5);
end


% set random number seed - comment out if desired.
%rand('state',9149);
%randn('state',9149);


for trials=1:ntrials, trials                %Loop through a number of trials
    for auxc = 1:length(coupling)      %The range of couplings
        for auxs = 1:length(delays)     %The range of lags to be scanned.
            rstate=rng;
            
            c=coupling(auxc);     % ********* COUPLING ***********
            scanlag=delays(auxs);  % ********* DELAYS ***********
            starttime = 0;
            myrand = rand(N,1); %random seed used to generate different IC/histories (variable is passed to @nrmlmass_hist)
            myddehandle = @nrlmass_dde_3;     %dde file handle
            
            
            if scanlag == 0
                mylags = [];                   %list of lag times required to evaluate the coupled dde
            else
                mylags = scanlag;
            end
            
            
            myhisthandle = @nrlmass_hist;    %handle for function that gives dde history
            myoptions = ddeset('RelTol', 1e-6, 'AbsTol', 1e-6);   %options structure, for adjusting e.g. error tolerance
            
            
            fprintf('Solving DDE for lag = %g, coupling = %g ... ', scanlag, c)
            
            if auxc==1
                %1. Run initial conditions
                if rand_ics
                    sol = ddesolver(myddehandle, mylags, myhisthandle, [starttime ictime], myoptions);
                else
                    sol = ddesolver(myddehandle, mylags, ics_sol, [starttime ictime], myoptions);  %solve the DDE
                end
            else
                sol=solrestart;
                sol.x=sol.x-sol.x(end)+ictime; % realign time
                sol.discont=[]; % ok??
            end
            
            starttime=ictime; endtime=starttime+ictime;
            
            mysavename = sprintf('%s_d%.8fms_%s_coupling%.4f_trial%i.mat', basename, scanlag, label, c,trials);  %filename for saving mat-files and figures
            
            
            %perturb=trials-1;   %add a perturbation?
            perturb=0;   %add a perturbation?
            rand_node=1;  % for saving in the case of no perturb
            rand_time=0;  % for saving in the case of no perturb (and for rand_node case)
            
            if perturb
                if option_node   %random node
                    rand_node=randperm(length(connectivity));
                    rand_node_i=rand_node(1:nds)*3-2;
                else              %non-random node but further initial integration of random length
                    rand_node_i=(non_random_nde*3)-2;
                    rand_time=ceil(rand*rand_length)+offset;
                    starttime=ictime; endtime=starttime+rand_time;
                    tic
                    sol=ddesolver(myddehandle, mylags, sol, [starttime endtime], myoptions);
                    toc
                    starttime=endtime; endtime=starttime+ictime;
                end
                sol.y(rand_node_i,end)=sol.y(rand_node_i,end).*perturbed;
            end
            
            checkpoint=0;
            % solve DDEs
            for segment=1:segments
                fprintf('d=%g,c=%g: seg %d of %d: ',scanlag,c,segment,segments);
                tic
                sol=ddesolver(myddehandle, mylags, sol, [starttime endtime], myoptions);
                starttime=endtime; endtime=starttime+ictime;
                toc
                
                % save checkpoint
                if ~mod(segment,50)
                    checkpoint=checkpoint+1;
                    fprintf('saving checkpoint %d\n',checkpoint)
                    % downsample
                    
                    %time = linspace(sol.x(1),sol.x(end),tfactor*(sol.x(end)-sol.x(1)));  %time grid on which the DE solution is evaluated            
                    time = sol.x(1):outdt:sol.x(end);
                    soln = deval(sol, time, 1:3:N*3-2); %the DE solution; keep only V

                    solrestart=sol;
                    keep=sol.x>=(sol.x(end)-ictime);
                    solrestart.x=solrestart.x(keep);
                    solrestart.y=solrestart.y(:,keep);
                    solrestart.yp=solrestart.yp(:,keep);

                    checkptstr=sprintf('_checkpoint%d.mat',checkpoint);

                    %save output
                    save([outdir filesep mysavename checkptstr],'solrestart','time','soln','rand_node','perturb','rand_time','normW','perturbed','ictime','segments','segment','outdt','checkpoint','coupling','delays');

                    % truncate because otherwise dde23 is too slow
                    sol=solrestart;
                end
                
            end      %segments
            
                
            % save last chunk

            % downsample
            %time = linspace(sol.x(1),sol.x(end),tfactor*(sol.x(end)-sol.x(1)));  %time grid on which the DE solution is evaluated            
            time = sol.x(1):outdt:sol.x(end);
            soln = deval(sol, time, 1:3:N*3-2); %the DE solution; keep only V

            solrestart=sol;
            keep=sol.x>(sol.x(end)-ictime);
            solrestart.x=solrestart.x(keep);
            solrestart.y=solrestart.y(:,keep);
            solrestart.yp=solrestart.yp(:,keep);

            checkptstr=sprintf('_checkpointEND.mat');
            checkpoint=checkpoint+1;

            %save output
            save([outdir filesep mysavename checkptstr],'solrestart','time','soln','rand_node','perturb','rand_time','normW','label','perturbed','ictime','segments','segment','outdt','checkpoint','coupling','delays','rstate','myrand');
            
        end         %lags
    end             %coupling
end                 %trials

% put the checkpoints back together
brainwaves_long_checkpointmerge(outdir)