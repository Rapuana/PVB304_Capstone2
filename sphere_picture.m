function sphere_picture(yp, t, cmap)
% Draw a picture of the brain at a time point 't', cmap refers to the
% colour map used, and is by default jet
if nargin<3
    cmap='jet';
end
figure
crange=quantile(yp(:),[0.001 0.999]); % may be needlessly slow for big arrays
% figure(555), clf don't need a figure for this one
colormap(cmap)
ax=axes;
caxis([-0.25 0])
colorbar

nnodes=size(yp,2);

loc=getloc(nnodes);

sphereradius = 0.5;

sp=add_sphere_size_internal(loc(:,1),loc(:,2),loc(:,3),sphereradius*ones(nnodes,1),yp(1,:));
drawnow
cdata=get(sp,'cdata');

for k=1:nnodes 
    cdata(:,(k-1)*22 + (1:21))=yp(t,k);
end
set(sp,'cdata',cdata)

drawnow
title(sprintf('Time = %.0f ms', t))

az = 0; el = 0; roll = -90;


view(az,el)
camroll(roll)

end



function hout=add_sphere_size_internal( x, y, z, s, c)
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
% James Roberts, QIMR Berghofer, 2014-2016

hold on;

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
h=surf(Xall,Yall,Zall,Call,'EdgeColor','none');
%light('Position',[15 -20 10])
daspect([1 1 1])
if nargout>0
    hout=h;
end
end