function brainwaves_vid_spheres(outfilename,yvid,fps, cmap, az, el, roll, time)

if exist(outfilename,'file')
    overwrite=input([strrep(outfilename,'\','\\') ' exists, overwrite? y/[n]: '],'s');
    if ~strcmp(overwrite,'y')
        return
    end
end

% vid content
% given as input, but need node locs
nnodes=size(yvid,2);

loc = getloc(nnodes);

if nnodes==112
    sphereradius=0.5; % mouse brain is small
else
    sphereradius=4;
end


% setup output video
writerObj = VideoWriter(outfilename,'MPEG-4');
writerObj.FrameRate=fps;
writerObj.Quality=100; % 1=min, 100=best
open(writerObj);

outwidth=1024;
outheight=576; 

hf=figure(888); clf
figpos=[50 50 outwidth outheight];
set(hf,'position',figpos,'color',[1 1 1],'paperPositionMode','auto'); 
ax=axes('visible','off','units','normalized','position',[0 0 1 1],'color','w');
view(az,el)
camroll(roll)

% Spheres
colormap(cmap)

if nargin < 8
    time = [500, length(yvid)];
end
for jj=time(1):time(2)
    if ~mod(jj,round(fps)*60)
        fprintf(1,'%d ',jj);
        if ~mod(jj,20*round(fps)*60)
            fprintf(1,'\n');
        end
    end
    %cla
    
    % spheres
    if jj==time(1)
%         crange=quantile(yvid(:),[0.001 0.999]);
        crange = [-0.25 0]
        caxis(crange)
        cbar=colorbar;
        ylabel(cbar,'V')
        sp=add_sphere_size_internal(ax,loc(:,1),loc(:,2),loc(:,3),sphereradius*ones(nnodes,1),yvid(1,:));
        cdata=get(sp,'cdata');
        th=text(0.0,0.95,sprintf('t = %g ms',jj-1),'units','normalized','fontsize',18,'horizontalalignment','left');
    else
    % % faster version (but maybe drawnow is rate-determining step on server)
        for k=1:nnodes
            cdata(:,(k-1)*22 + (1:21))=yvid(jj,k);
        end
        set(sp,'cdata',cdata)
        set(th,'string',sprintf('t = %g ms',jj-1))
    end
    
    %set(ax,'cameraviewangle',5)
    %set(ax,'visible','off');
    
    % write output frame
    writeVideo(writerObj,getframe(hf));
    %img = hardcopy(h, '-dzbuffer', '-r0');
    %writeVideo(writerObj, im2frame(img));
end
close(writerObj);

end



function hout=add_sphere_size_internal(ax, x, y, z, s, c)
% input:    x   Centrepoint of the sphere in the x direction
%           y   Centrepoint of the sphere in the y direction
%           z   Centrepoint of the sphere in the z direction
%           s   Size for each sphere
%           c   Value by which to color the sphere
%
% output: (optional) handles to the surfs
% 
%
% Anton Lord, UQ 2010
% James Roberts, QIMR Berghofer, 2014-2017

hold(ax,'on');

n = 20; % number of faces each sphere has
c = double(c);
[x0,y0,z0] = sphere(n);

nspheres=length(x);
Xall=nan(n+1,nspheres*(n+2)); % each sphere is an (n+1)-by-(n+1) matrix
Yall=Xall;
Zall=Xall;
Call=Xall;

for j = 1:length(x)
    if size(c,1) == 1
        intensity = zeros(n+1)+c(j); 
    end
    Xall(:,(j-1)*(n+2)+(1:(n+1)))=x0*s(j)+x(j);
    Yall(:,(j-1)*(n+2)+(1:(n+1)))=y0*s(j)+y(j);
    Zall(:,(j-1)*(n+2)+(1:(n+1)))=z0*s(j)+z(j);
    Call(:,(j-1)*(n+2)+(1:(n+1)))=intensity;
end
h=surf(ax,Xall,Yall,Zall,Call,'EdgeColor','none');
%light('Position',[15 -20 10])
daspect(ax,[1 1 1])
if nargout>0
    hout=h;
end
end
