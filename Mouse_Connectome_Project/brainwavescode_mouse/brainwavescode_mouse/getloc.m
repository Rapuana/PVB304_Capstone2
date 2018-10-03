function loc=getloc(nnodes)
%  loc = getloc(nnodes)
% get nnodes-by-3 node coordinate matrix given the number of nodes nnodes
% in the parcellation

switch nnodes
    case 513
        if isunix % if on our cluster with ldrive mounted in ~/ldrive
            %load('~/ldrive/Lab_MichaelB/Rich/connectomes/distmat','loc')
            loc=load('~/ldrive/Lab_MichaelB/jamesR/networks/513COG.mat'); loc=loc.COG;
        else
            %load('L:\Lab_MichaelB\Rich\connectomes\distmat','loc')
            loc=load('L:\Lab_MichaelB\jamesR\networks\513COG.mat'); loc=loc.COG;
        end
    case 512
        if isunix % if on our cluster with ldrive mounted in ~/ldrive
            loc=load('~/ldrive/Lab_MichaelB/alistair/newparcinfo/COGnew.mat'); loc=loc.COG;
        else
            loc=load('L:\Lab_MichaelB\alistair\newparcinfo\COGnew.mat'); loc=loc.COG;
        end
    case 58
        loc=load('MyAtlas_n58.mat'); loc=loc.MyAtlas.Centroids*1000; % into mm
    case 84
        loc=load('fs_default_xyz'); loc=loc.fs_default_xyz;
    case 82
        loc=load('fs_default_xyz'); loc=loc.fs_default_xyz;
        loc=loc([1:34,36:83],:); % delete cerebellum
    case 90
        loc=load('~/ldrive/Lab_MichaelB/jamesR/networks/fernando/connectomes/ferHCP_aal.mat');
        loc=loc.ferHCP.loc;
    case 112
        %loc=load('L:\Lab_MichaelB\jamesR\networks\mouse\rubinovmouse.mat');
		loc=load('rubinovmouse.mat');
        loc=loc.rubinovmouse.XYZ/1000; % into mm?
    case 998
        if isunix % if on our cluster with ldrive mounted in ~/ldrive
            loc=load('~/ldrive/Lab_MichaelB/jamesR/networks/Hagmann_PLoSBiol_2008_group_mean_region_xyz_centers.txt');
        else
            loc=load('L:\Lab_MichaelB\jamesR\networks\Hagmann_PLoSBiol_2008_group_mean_region_xyz_centers.txt');
        end
end