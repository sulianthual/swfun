function animatesw(filename,frameskip,dolevels,x,y,t,u,v,p)
% 
% This function creates an time-animation of maps x-y of fields {u,v,p}(x,y,t)
% It contours the animation and saves it in a .gif file
%
% Input parameters:
%  - filename: name of .gif file, for example filename='shallowwater';
%  - frameskip: how often to draw, for example frameskip=10 will draw only every 10 timesteps (take frameskip=1 firsT)
%  - dolevels: =1 keep same contours levels each frame or =0 dont (try dolevels=1 first). For dolevels=0 there is no colorbars.
%  - x,y,t: the axis in array form
%  - u,v,p: the fields in multidimensional arrays form
%  
% (by S.Thual, Fudan 2019, class "GFD, theory and modeling")
%
% Ensure x,y,t in row or array are treated the same
Nx=length(x);
Ny=length(y);
Nt=length(t);
xx=zeros(1,Nx); xx(1,:)=x(:);
yy=zeros(1,Ny); yy(1,:)=y(:);
tt=zeros(1,Nt); tt(1,:)=t(:);

% levels for contours
if dolevels==1
umin=min(u(:)); umax=max(u(:)); ulevs=(0:0.1:1)*(umax-umin)+umin;
vmin=min(v(:)); vmax=max(v(:)); vlevs=(0:0.1:1)*(vmax-vmin)+vmin;
pmin=min(p(:)); pmax=max(p(:)); plevs=(0:0.1:1)*(pmax-pmin)+pmin;
end

% Contour with time loop
figure(2);
for k=1:frameskip:Nt-1
clf; 
%
subplot(1,3,1)
if dolevels==1
contourf(xx,yy,squeeze(u(:,:,k))',ulevs,'linestyle','none'); caxis([umin,umax])
else
contourf(xx,yy,squeeze(u(:,:,k))','linestyle','none');
end
title('u'); xlabel('x'); ylabel('y'); 
if dolevels==1; colorbar('location','southoutside'); 
else colorbar('location','southoutside','visible','off'); end
%
subplot(1,3,2)
if dolevels==1
contourf(xx,yy,squeeze(v(:,:,k))',vlevs,'linestyle','none'); caxis([vmin,vmax])
else
contourf(xx,yy,squeeze(v(:,:,k))','linestyle','none');
end
title('v'); xlabel('x'); ylabel('y');
if dolevels==1; colorbar('location','southoutside'); 
else colorbar('location','southoutside','visible','off'); end

%
subplot(1,3,3)
if dolevels==1
contourf(xx,yy,squeeze(p(:,:,k))',plevs,'linestyle','none'); caxis([pmin,pmax])
else
contourf(xx,yy,squeeze(p(:,:,k))','linestyle','none');
end
title(strcat(['p (t=',num2str(tt(k)),')'])); xlabel('x'); ylabel('y');
if dolevels==1; colorbar('location','southoutside'); 
else colorbar('location','southoutside','visible','off'); end
%
% Save to gif file
drawnow
frame = getframe(2); im = frame2im(frame); [imind,cm] = rgb2ind(im,256);
if k==1; 
imwrite(imind,cm,strcat([filename,'.gif']),'gif', 'Loopcount',inf);
else; 
imwrite(imind,cm,strcat([filename,'.gif']),'gif','WriteMode','append'); 
end

end% end time loop

