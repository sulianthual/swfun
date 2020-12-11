% By Sulian Thual
%
% Solve the shallow-water equations:
% dtu-fv=-dxp
% dtv+fu=-dyp
% dtp+G(dxu+dyv)=0
%
% vertical mode n: un=u,vn=v,pn=p/rho,cn=sqrt(G)
% one layer: u1=u,v1=v,h1=p/g, H=G/g
% reduced gravity: u2=u, v2=v, h=p/g*, H2=G/g*, g*=(rho2-rho1)*g/rho1
% 
% Finite difference methods: 
% -explicit Euler scheme in time
% -centered scheme in space
%
% test for different initial/boundary conditions

clear all;

disp('1) initialization');disp(' ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Parameters

c=0.1; % phase speed (km.day-1)
G=c^2;% always positive 


bdy_x=1; % bdy conditions east/west (=1 perio, =2 reflecting)
bdy_y=1; % bdy conditions north/south (=1 perio, =2 reflecting)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Case of Grid/flow for evaluation

% All midterm examples have same:
% initial condition (u=0,v=0,p=centered gaussian)
% periodic conditions in x
% reflecting conditions at ymin, ymax
% c=0.1 (G=0.01)

% midterm Task 1: periodic, f=0
% cflx,cfly=0.02,0.02
% f=0
if 1==1
filename = 'midterm1'; 
Lx=500; Nx=200; dx=2.5;
Ly=500; Ny=200; dy=2.5;
Lt=5000; Nts=100; Ntskip=100;  dt=0.5; Nt=Ntskip*Nts;
frameskip=1; % for animatesw
 % save Nts points each Ntskip timesteps
f=zeros(1,Ny); % Coriolis parameter
cfl=c*dt/dx; % CFL
aax=0.002; xxa=Lx/2; aay=0.002; yya=Ly/2;
bdy_x=1; bdy_y=1;% bdy conditions
end

% midterm but coarser
if 0==1
filename = 'midterm1a'; 
Nres=10; % resolution factor (>1 for coarser, <1 for finer, always integer or 1/)
Lx=500; Nx=200/Nres; dx=Nres*2.5;
Ly=500; Ny=200/Nres; dy=Nres*2.5;
Lt=5000; Nts=100; Ntskip=100/Nres; dt=Nres*0.5;  Nt=Ntskip*Nts;
frameskip=1; % for animatesw
f=zeros(1,Ny); % Coriolis parameter
cfl=c*dt/dx; % CFL
aax=0.002; xxa=Lx/2; aay=0.002; yya=Ly/2;
bdy_x=1; bdy_y=1;% bdy conditions
end

% midterm but save all timesteps
if 0==1%
filename = 'midterm1b'; 
Lx=500; Nx=200; dx=2.5;
Ly=500; Ny=200; dy=2.5;
Lt=5000; Nts=10000; Ntskip=1;  dt=0.5; Nt=Ntskip*Nts;
frameskip=100; % for animatesw
 % save Nts points each Ntskip timesteps
f=zeros(1,Ny); % Coriolis parameter
cfl=c*dt/dx; % CFL
aax=0.002; xxa=Lx/2; aay=0.002; yya=Ly/2;
bdy_x=1; bdy_y=1;% bdy conditions
end

% test f=cte
if 0==1
filename = 'midterm2'; 
Lx=500; Nx=200; dx=2.5;
Ly=500; Ny=200; dy=2.5;
Lt=5000; Nts=100; Ntskip=100; dt=0.5;  Nt=Ntskip*Nts;
frameskip=1; % for animatesw
f=zeros(1,Ny)+10e-2; 
cfl=c*dt/dx; % CFL
aax=0.002; xxa=Lx/2; aay=0.002; yya=Ly/2;
bdy_x=1; bdy_y=1;% bdy conditions
end

% test f=b(y-yo)
if 0==1
filename = 'midterm3'; 
Lx=500; Nx=200; dx=2.5;
Ly=500; Ny=200; dy=2.5;
Lt=5000; Nts=100; Ntskip=100; dt=0.5; Nt=Ntskip*Nts;
frameskip=1; % for animatesw
f=((0:1:Ny-1)*dy-Ly/2)*2*10e-5;
cfl=c*dt/dx; % CFL
aax=0.002; xxa=Lx/2; aay=0.002; yya=Ly/2;
bdy_x=1; bdy_y=1;% bdy conditions
end

% test f=cos(...)
if 0==1
filename = 'midterm4'; 
Lx=500; Nx=200; dx=2.5;
Ly=500; Ny=200; dy=2.5;
Lt=5000;  Nts=100; Ntskip=100; dt=0.5;  Nt=Ntskip*Nts;
frameskip=1; % for animatesw
f=-2*(10e-2)*cos((0:1:Ny-1)*dy*pi/Ly);
cfl=c*dt/dx; % CFL
aax=0.002; xxa=Lx/2; aay=0.002; yya=Ly/2;
bdy_x=1; bdy_y=1;% bdy conditions
end

% Fast for tests
if 0==1
filename = 'fast'; 
Lx=500; Nx=200; dx=2.5;
Ly=500; Ny=200; dy=2.5;
Lt=5000; Nts=30; Ntskip=200;  dt=0.5; Nt=Ntskip*Nts;
frameskip=1; % for animatesw
 % save Nts points each Ntskip timesteps
f=-2*(10e-2)*cos((0:1:Ny-1)*dy*pi/Ly); 
cfl=c*dt/dx; % CFL
aax=0.002; xxa=Lx/2; aay=0.002; yya=Ly/2;
bdy_x=1; bdy_y=1;% bdy conditions
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check on numerical stability (CFL,Nx*Ny*Nt)
if     1==1
disp(strcat(['Lx(km),Ly(km),Lt(days)=',num2str(Lx),',',num2str(Ly),',',num2str(Lt)]));
disp(strcat(['dx (km),dy (km),dt (days)=',num2str(dx),',',num2str(dy),',',num2str(dt)]));
disp(strcat(['Nx,Ny,Nt=',num2str(Nx),',',num2str(Ny),',',num2str(Nt)]));
disp(strcat(['CFLx,CFLy=',num2str(c*dt/dx),',',num2str(c*dt/dy)]));
disp(' ');
%  stop
end
%
x=(0:1:Nx-1)*dx; % grid
y=(0:1:Ny-1)*dy; % grid
t=(0:1:Nts-1)*dt*Ntskip; % grid

us=zeros(Nx,Ny,Nts); % array u
vs=zeros(Nx,Ny,Nts); % array v
ps=zeros(Nx,Ny,Nts); % array p

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Parameters again

% Initial conditions
u=zeros(Nx,Ny)+5;
v=zeros(Nx,Ny)+5;
p=exp(-aax*(x-xxa).^2)'*exp(-aay*(y-yya).^2); 
ue=u*0;
ve=v*0;
pe=p*0; 

if 0==1
figure(2); clf;
subplot(1,3,1); contourf(x,y,squeeze(us(:,:,1))','linestyle','none'); 
title('initial conditions: u'); xlabel('x'); ylabel('y');
subplot(1,3,2); contourf(x,y,squeeze(vs(:,:,1))','linestyle','none'); 
title('initial conditions: v'); xlabel('x'); ylabel('y');
subplot(1,3,3); contourf(x,y,squeeze(ps(:,:,1))','linestyle','none'); 
title('initial conditions: p'); xlabel('x'); ylabel('y');
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Technical

disp('2) Time integration'); disp(' ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Integrate the model

%  1) starting
% dtu-fv=-dxp
% dtv+fu=-dyp
% dtp+G(dxu+dyv)=0
% 2)  df/dx|i=(fi+1-fi-1)/2/dx centered
% dtu|i=fv-dxp=fvi -(pi+1-pi-1)/2/dx
% dtv|j=-fu-dyp=-fuj -(pj+1-pj-1)/2/dy
% dtp|i,j=-G(dxu+dyv)=-G*(ui+1-ui-1)/2/dx-G*(vj+1-vj-1)/2/dy
% 3) fk+1=fk+ dt*Qk: Euler explicit for Q=dtf
%  ui,j,k+1=ui,j,k +dt*f(j)*vi,j,k -dt*(pi+1,j,k - pi-1,j,k)/2/dx
%  vi,j,k+1=vi,j,k -dt*f(j)*ui,j,k -dt*(pi,j+1,k - pi,j-1,k)/2/dy
%  pi,j,k+1=pi,j,k -dt*G*(ui+1,j,k - ui-1,j,k)/2/dx - dt*G*(vi,j+1,k - vi,j-1,k)/2/dy 

% Bdy condition
% - periodic bdy: obvious 
% - reflecting bdy y=cte: u=0,v=0,dyp=0 (first order)
% - open bdy x=cte: buffer layer with diffusion, increase viscosity progressively
%
% Time integration
for ks=1:Nts
%
disp(strcat([num2str(fix(ks/Nts*100)),'%']));
% save
us(:,:,ks)=u(:,:);
vs(:,:,ks)=v(:,:);
ps(:,:,ks)=p(:,:);
%
for k=1:Ntskip

% Solve within domain
for i=2:Nx-1
for j=2:Ny-1
ue(i,j)=u(i,j) +dt*f(j)*v(i,j) -dt*(p(i+1,j) - p(i-1,j))/2/dx;
ve(i,j)=v(i,j) -dt*f(j)*u(i,j) -dt*(p(i,j+1) - p(i,j-1))/2/dy;
pe(i,j)=p(i,j) -dt*G*(u(i+1,j) - u(i-1,j))/2/dx - dt*G*(v(i,j+1) - v(i,j-1))/2/dy;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edge Periodic in x
if bdy_x==1
% Solve west boundary i=1
for j=2:Ny-1
ue(1,j)=u(1,j) +dt*f(j)*v(1,j) -dt*(p(2,j) - p(Nx,j))/2/dx; % periodic x
ve(1,j)=v(1,j) -dt*f(j)*u(1,j) -dt*(p(2,j+1) - p(Nx,j-1))/2/dy;
pe(1,j)=p(1,j) -dt*G*(u(2,j) - u(Nx,j))/2/dx - dt*G*(v(1,j+1) - v(1,j-1))/2/dy;
end

% Solve east boundary i=Nx
for j=2:Ny-1
ue(Nx,j)=u(Nx,j) +dt*f(j)*v(Nx,j) -dt*(p(1,j) - p(Nx-1,j))/2/dx;% periodic x
ve(Nx,j)=v(Nx,j) -dt*f(j)*u(Nx,j) -dt*(p(Nx,j+1) - p(Nx,j-1))/2/dy;
pe(Nx,j)=p(Nx,j) -dt*G*(u(1,j) - u(Nx-1,j))/2/dx - dt*G*(v(Nx,j+1) - v(Nx,j-1))/2/dy;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edge Periodic in y
if bdy_y==1
% Solve south boundary j=1
for i=2:Nx-1
ue(i,1)=u(i,1) +dt*f(1)*v(i,1) -dt*(p(i+1,1) - p(i-1,1))/2/dx;% periodic y
ve(i,1)=v(i,1) -dt*f(1)*u(i,1) -dt*(p(i,2) - p(i,Ny))/2/dy;
pe(i,1)=p(i,1) -dt*G*(u(i+1,1) - u(i-1,1))/2/dx - dt*G*(v(i,2) - v(i,Ny))/2/dy;
end

% Solve north boundary j=Ny
for i=2:Nx-1
ue(i,Ny)=u(i,Ny) +dt*f(Ny)*v(i,Ny) -dt*(p(i+1,Ny) - p(i-1,Ny))/2/dx;% periodic y
ve(i,Ny)=v(i,Ny) -dt*f(Ny)*u(i,Ny) -dt*(p(i,1) - p(i,Ny-1))/2/dy;
pe(i,Ny)=p(i,Ny) -dt*G*(u(i+1,Ny) - u(i-1,Ny))/2/dx - dt*G*(v(i,1) - v(i,Ny-1))/2/dy;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Corners Periodic in x, Periodic in y
if bdy_x==1
if bdy_y==1
% Solve corner i=1,j=1
ue(1,1)=u(1,1) +dt*f(1)*v(1,1) -dt*(p(2,1) - p(Nx,1))/2/dx;% periodic x|periodic y
ve(1,1)=v(1,1) -dt*f(1)*u(1,1) -dt*(p(1,2) - p(i,Ny))/2/dy;
pe(1,1)=p(1,1) -dt*G*(u(2,1) - u(Nx,1))/2/dx - dt*G*(v(1,2) - v(1,Ny))/2/dy;

% Solve corner i=Nx,j=1
ue(Nx,1)=u(Nx,1) +dt*f(1)*v(Nx,1) -dt*(p(1,1) - p(Nx-1,1))/2/dx;% periodic x|periodic y
ve(Nx,1)=v(Nx,1) -dt*f(1)*u(Nx,1) -dt*(p(Nx,2) - p(Nx,Ny))/2/dy;
pe(Nx,1)=p(Nx,1) -dt*G*(u(1,1) - u(Nx-1,1))/2/dx - dt*G*(v(Nx,2) - v(Nx,Ny))/2/dy;

% Solve corner i=1,j=Ny
ue(1,Ny)=u(1,Ny) +dt*f(Ny)*v(1,Ny) -dt*(p(2,Ny) - p(Nx,Ny))/2/dx;
ve(1,Ny)=v(1,Ny) -dt*f(Ny)*u(1,Ny) -dt*(p(1,1) - p(1,Ny-1))/2/dy;
pe(1,Ny)=p(1,Ny) -dt*G*(u(2,Ny) - u(Nx,Ny))/2/dx - dt*G*(v(1,1) - v(1,Ny-1))/2/dy;

% Solve corner i=Nx,j=Ny
ue(Nx,Ny)=u(Nx,Ny) +dt*f(Ny)*v(Nx,Ny) -dt*(p(1,Ny) - p(Nx-1,Ny))/2/dx;
ve(Nx,Ny)=v(Nx,Ny) -dt*f(Ny)*u(Nx,Ny) -dt*(p(Nx,1) - p(Nx,Ny-1))/2/dy;
pe(Nx,Ny)=p(Nx,Ny) -dt*G*(u(1,Ny) - u(Nx-1,Ny))/2/dx - dt*G*(v(Nx,1) - v(Nx,Ny-1))/2/dy;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edges Reflecting in x
if bdy_x==2
% Solve west boundary i=1
for j=1:Ny 
ue(1,j)=0; % wall in x
ve(1,j)=0;
pe(1,j)=p(2,j);
end

% Solve east boundary i=Nx
for j=1:Ny 
ue(Nx,j)=0;% wall in x
ve(Nx,j)=0;
pe(Nx,j)=p(Nx-1,j);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edges Reflecting in y
if bdy_y==2
% Solve south boundary j=1
for i=1:Nx 
ue(i,1)=0;% wall in y
ve(i,1)=0;
pe(i,1)=p(i,2);
end

% Solve north boundary j=Ny
for i=1:Nx 
ue(i,Ny)=0;% wall in y
ve(i,Ny)=0;
pe(i,Ny)=p(i,Ny-1);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edges Open in y (wave equation normal to bdy)
% dtf+c*dyf=0, c towards exterior, for f=u,v,p
if 0==1

% Solve south boundary j=1
for i=1:Nx 
ue(i,1)=u(i,1) +sqrt(G)*dt/dx*(u(i,1) -u(i,2));% open in y (wave eq)
ve(i,1)=v(i,1) +sqrt(G)*dt/dx*(v(i,1) -v(i,2));
pe(i,1)=p(i,1) +sqrt(G)*dt/dx*(p(i,1) -p(i,2));
end

% Solve north boundary j=Ny
for i=1:Nx
ue(i,Nx)=u(i,Nx) -sqrt(G)*dt/dx*(u(i,Nx) -u(i,Nx-1));% open in y (wave eq)
ve(i,Nx)=v(i,Nx) -sqrt(G)*dt/dx*(v(i,Nx) -v(i,Nx-1));
pe(i,Nx)=p(i,Nx) -sqrt(G)*dt/dx*(p(i,Nx) -p(i,Nx-1));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update
u(:,:)=ue;
v(:,:)=ve;
p(:,:)=pe;

end% k

end% ks

disp('3) Analysis'); disp(' ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Show results

% Animate and save in self-contained function
animatesw(filename,frameskip,0,x,y,t,us,vs,ps);


% Hovmollers at y=Ly/2
if 0==1
figure(4); clf; 
pass=fix(Ny/2); 
subplot(1,3,1)
contourf(x,t,squeeze(us(:,pass,:))','linestyle','none'); 
title('u(x,Ly/2,t)'); xlabel('x'); ylabel('t');
subplot(1,3,2)
contourf(x,t,squeeze(vs(:,pass,:))','linestyle','none'); 
title('v(x,Ly/2,t)'); xlabel('x'); ylabel('t');
subplot(1,3,3)
contourf(x,t,squeeze(ps(:,pass,:))','linestyle','none'); 
title('p(x,Ly/2,t)'); xlabel('x'); ylabel('t');
end

% gifview shortcuts:
%  a:Toggle between animation and slideshow mode.
%  Space: Go to the next frame.
%  b: Go to the previous frame.
%  r: Go to the first frame.
%  >: Go to the last frame.
%  u:Toggle between normal and unoptimized mode.
%  q: Quit gifview.
