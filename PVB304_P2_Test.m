%% Testing for PVB304 - Project 2




t_perturb = 1304;

segments = 50;

ictime = 50;





t_perturb/ictime

26 * ictime

% Segment in which the perturbation occurs
perturb_seg = floor(t_perturb/ictime);
;
for segment = 1:segments
    if segment == perturb_seg
        % This endtime is for when there is no perturbation
        endtime = t_perturb - perturb_seg * ictime - 1;
        
        % DDE23 solver here
        
        starttime = endtime + 1;
        endtime = starttime + t_length;
        Perturbed = strength_perturb; % Hopefully this works. 
        % Essentially I am just adding this value onto the ODE over the
        % entirety of the solution. If not I have a few other ideas. 
        
        % DDE23 solver here (PERTURBED)
        
    end
    
end
        