function mysavename = sam_perturbations_BW(normW, coupling, delays, ...
    label, outdir, t_perturb, strength_perturb, node)
% A function adapted from James Roberts and Anton Lord to apply
% perturbations to the mouse connectome. This function will be written with
% the mouse connectome specifically in mind and thus may be inefficient on
% larger connectomes. 
%
% Author: Samuel Dudley ~ 
% Queensland Institute of Medical Research (Visiting Student)
% Queensland University of Technology (BSc. & BMath.)
%
% Adapted from: James Roberts,
%               Anton Lord.
%
% sam_perturbations_BW
% normW: Connectivity Matrix
% coupling: normalised coupling strength from 0 to 1
% delays: time of delay between solves, (in ms)
% label: string -> Code label. example: 'mouse'
% outdir: string -> Directory to save files
% perturb: true or false value 
% t_perturb: time at which perturbation begins
% t_length: time for which perturbation occurs REVISION -- REMOVED
% strength_perturb: Strength of the perturbation
% node: node in which we perturb

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% NOTES % % NOTES % % NOTES % % NOTES % % NOTES % % NOTES % % NOTES %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DDE23 is a Delayed Differential Equation solver using Runge Kutta's 2nd
% and 3rd order methods. Higher order methods are very computationally
% heavy and not required.

%%% sol = DDE23(ddefun,lags,history,tspan,options)
% ddefun DDE Function which has the form
% dydt = ddefun(t,y,Z)
% whereby t is the current time, 
% y is a column vector approximating y(t) and Z(:,j) approximates
% y(t-tau_j) for delay tau_j = lags(j). The output is a column vector
% corresponding to f(t,y(t), y(t-tau_1),...,y(t-tau_k)).
% lags DELAYS which has length(lags) = k and is a vector of constant
% positive delays. 
% tspan refers to the interval of integration from t_0=tspan(1) to
% t_f=tspan(end) such that t_i < t_i+1
% sol SOLUTION 
% the structure of the solution given by dde23. 
% sol.x is the Mesh selected by dde23 to solve over
% sol.y is the approximation to y(x) at points in sol.x
% sol.yp is the approximation to y'(x) at points in sol.x
% sol.solver = dde23 (the solver name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% SOLVER % % SOLVER % % SOLVER % % SOLVER % % SOLVER % % SOLVER %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set global variables which we share with DDE solver

global V1 V2 V3 V4 V5 V6 V7 gCa gK gL VK VL VCa I b ani aei aie aee phi ...
    V8 V9 gNa VNa ane nse ...
    rnmda N CM vs c k_in numddevars myrand

numddevars = 3;


%%% Set up values
% Set the timings for the solution. 
ictime = 50;        % initial conditions for the integration
segments = 1;      % the number of segments we split the integration into
connectivity = normW;


% Logic Tests
if nargin > 5 % We check to see if we are running perturbations.
    % If so check the perturbation code
    if t_perturb >= ictime*segments 
        error(['Perturbing at time outside range, ' ...
               't_perturb must be less than ictime*segments'])
    elseif t_perturb < ictime
        error(['Perturbing at a time too early, ' ...
               't_perturb must be greater than ictime'])
    end
    if node > length(normW)
        error('You are trying to perturb a node that does not exist')
    end
    % Set the perturb logical
    perturb = true;
else
    % Set the perturb logical
    perturb = false;
end


% Set the output timestep (in ms)
outdt = 1;  % Default = 1

% Set a name for the output of the filename base
basename = ['Run_ictime' num2str(ictime) '_seg' num2str(segments)];


% Set directory 
mkdir(outdir);



% Set up Sparse Matrix

% Due to the size of the mouse connectome this may not be worth it for the
% mouse connectome from a computational stand point. Regardless we shall
% leave it.
% Initialise sparse matrix 
CM = sparse(connectivity);

% Take size of connectivity matrix
[N,~] = size(connectivity); % Should be square.

% Sum all the columns of the connectivity matrix and transpose into col vec
k_in = sum(CM)';

% Set up model Parameters

V1 = -0.01; V2 = 0.15; V3 = 0   ; V4 = 0.3; V5 = 0; 
V6 =  0.65; V7 = 0   ; V8 = 0.15; V9 = 0.3;
gCa = 1; gK = 2.0; gL = 0.5; gNa = 6.7;
VK = -0.7; VL = -0.5; I = 0.3; b = 0.1; phi = 0.7; VNa = 0.53; VCa = 1;
ani = 0.4; vs = 1; aei = 2; aie = 2; aee = 0.36; ane = 1; rnmda = 0.25;


% Noise, coupling, modulation
% NOTE %
% Our model does not incorporate noise but regardless we shall leave this
% global variable here for future implementations.
nse = 0;




%%% Now we are onto the actual solving

% We have removed loops for 'trials', 'coupling' and 'delay'. 

rstate = rng;           % Do we need this?
c = coupling;           % Coupling GP
scanlag = delays;       % Delays

starttime = 0;          % Time we start solve at
myrand = rand(N,1);     % random seed used to generate different IC 
                        % histories  which is parsed to $nrlmass_hist
myddehandle = @nrlmass_dde_3; % DDE file handle

% Set the solver to dde23.
ddesolver=@(varargin) dde23(varargin{:});


if scanlag == 0
    mylags = []; 
else
    mylags = scanlag;
end

% Set handle for function that looks after dde history
% History from before t_0
myhisthandle = @nrlmass_hist; 
myoptions = ddeset('RelTol', 1e-6, 'AbsTol', 1e-6); 

fprintf('Solving DDE for lag = %g, coupling = %g ... ', scanlag, c)

% Assume random initial conditions
sol = ddesolver(myddehandle, mylags, myhisthandle, ...
                [starttime ictime], myoptions);

% Run over the initial conditions. 
starttime = ictime; endtime = starttime + ictime; 

% Set the save name for the mat-files and figures
mysavename = ...
    sprintf('%s_d%.3fms_%s_coupling%.3f_Pert%g_Pert_time%g.mat', ....
                basename, scanlag, label, c, perturb, t_perturb);

%%% This is where we add the perturbations

% Arguments used
% t_length <- DELETED
% t_perturb
% strength_perturb
% node


% Perturb arg is set based on nargin

% We need to determine which segment we run the perturbation in.
if perturb % This will run if perturb is active
    % Determine segment in which perturb will occur. 
    perturb_seg = floor(t_perturb/ictime);
end

for segment = 1:segments
    if segment == perturb_seg % Check to see if we add perturbation in this
        % loop. If so we change the endtime
        % STARTTIME IS AS USUAL
        endtime = t_perturb; % Set endtime to the time in which perturb 
        % occurs. 
        
        % Let's have a sanity check on our coding here. 
        if endtime - starttime > ictime
            warning('PERTURBATION CODE WARNING:')
            warning('The selected endtime is larger than expected...')
        end
        
        fprintf('PERTURBATION SEGMENT: d=%g, c=%g, segment %d of %d', ...
            scanlag, c', segment, segments)
        fprintf('Perturbing node %g, at time %g with a strength of %g', ...
            node, t_perturb, strength_perturb)
        fprintf('With shorter segment, ictime = %g,\n', ...
            diff([starttime endtime]))
        tic % Time the iteration
        % Run the DDE solver as usual, stopping at perturbation point.
        sol = ddesolver(myddehandle, mylags, sol, ...
                        [starttime endtime], myoptions);
                    
        % Now we add the perturbation to the end of the desired node.
        sol.y(node, end) = sol.y(node, end) + strength_perturb;
        
        % Now we set the start time for the next segment.
        starttime = endtime; % Set this differently 
        toc% Iteration timing
        
    else
        if segment == perturb_seg+1
            fprintf('Longer segment, ictime = %g', ...
                diff([starttime endtime]))
        end
        fprintf('d=%g, c=%g,  segment %d of %d\n', ...
            scanlag, c, segment, segments)
        tic % Time iterations
        sol = ddesolver(myddehandle, mylags, sol, ...
                        [starttime endtime], myoptions);
        starttime = endtime; endtime = starttime + ictime;
        toc % Iteration timing
    end
end

% Save values
time = sol.x(1):outdt:sol.x(end);
soln = deval(sol, time, 1:3:N*3-2);

solrestart = sol;
keep = sol.x > (sol.x(end) - ictime);
solrestart.x = solrestart.x(keep);
solrestart.y = solrestart.y(:,keep);
solrestart.yp = solrestart.yp(:,keep);


% SAVING
% Save stuff here.

save([outdir filesep mysavename],                                       ...
    'solrestart','time','soln','node','perturb','t_perturb',            ...
    'normW','label','strength_perturb','ictime','segments','segment',   ...
    'outdt','coupling','delays','rstate','myrand', 'perturb_seg');


end % END FUNCTION