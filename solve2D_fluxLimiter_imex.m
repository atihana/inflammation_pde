function solve2D_fluxLimiter_imex(weight,fname)

%m3chem2d_cil_f
%model 3 2D implementation
%all parameters dimensionless
%periodic boundary conditions
%vectorized form with only column matrices
%space discretization only --> solving ODEs system with ode15s
%jacobian matrix as input to improve computation time
%different initial conditions considering centred blob
%--> no use of damage function --> alpha parameter not used

%ALL diffusion constants are considered
%chemotaxis rates for n and m

%0.310848 0.009237 21.194060 0.262276 0.025078 1.281019

%clear all
%close all
clc
%profile on
nX=60;
nY=60;
Lx=linspace(0,1,nX);    %current code: N=60-->52 minutes                 
Ly=linspace(0,1,nY);

[x,y]=meshgrid(Lx,Ly);  %asse cartesiano
dx = Lx(2)-Lx(1);              %dx=dy
L=length(Lx);           %=N=number of nodes

dt=0.25;
tmax=1000;         
t=0:dt:tmax;
numTimesteps=numel(t);

BC='Periodic';

% Set parameter values
%s=setParameterValues('phi',phi,'nu',nu);
s=setParameterValues_ParamSet2('phi',0.1,'nu',0.075);
s.ni=s.nu;

% Find corresponding steady states
y0=[3.9406,0.0770,38.2178,0.3822,0.0294]';
ICs=fsolve(@(y)model3_odes(0,y,s,0),y0);
n0=weight*ICs(1);
a0=weight*ICs(2);
m0=weight*ICs(3);
c0=weight*ICs(4);
g0=weight*ICs(5);

% Set initial conditions
damageRadius=0.25;
noiseAmplitude=0;

n=zeros(L,L);
n(((x-0.5).^2+(y-0.5).^2)<damageRadius^2)=n0;
n=n(:);
n=n.*(1+noiseAmplitude*randn(nX*nY,1));

a=zeros(L,L);
a(((x-0.5).^2+(y-0.5).^2)<damageRadius^2)=a0;
a=a(:);
a=a.*(1+noiseAmplitude*randn(nX*nY,1));

m=zeros(L,L);
m(((x-0.5).^2+(y-0.5).^2)<damageRadius^2)=m0;
m=m(:);
m=m.*(1+noiseAmplitude*randn(nX*nY,1));

c=zeros(L,L);
c(((x-0.5).^2+(y-0.5).^2)<damageRadius^2)=c0;
c=c(:);
c=c.*(1+noiseAmplitude*randn(nX*nY,1));

g=zeros(L,L);              
g(((x-0.5).^2+(y-0.5).^2)<damageRadius^2)=g0;
g=g(:);
g=g.*(1+noiseAmplitude*randn(nX*nY,1));

s.dam=0;

%boundary conditions
%diffusion
s.dmat=sparse(getDiffMatrix(L,dx,BC));
%chemotaxis
switch BC
    case 'Periodic'
        for i=1:L^2
                if(mod(i,L)==1), left(i)=i-1+L; else left(i)=i-1; end
                if(mod(i,L)==0), right(i)=i+1-L; else right(i)=i+1; end
                if(i>(L-1)*L),   up(i)=i-(L-1)*L; else up(i)=i+L; end
                if(i<L+1),       down(i)=(L-1)*L+i; else down(i)=i-L; end
        end
    case 'Neumann'
        for i=1:L^2
                if(mod(i,L)==1), left(i)=i+1; else left(i)=i-1; end
                if(mod(i,L)==0), right(i)=i-1; else right(i)=i+1; end
                if(i>(L-1)*L),   up(i)=i-L; else up(i)=i+L; end
                if(i<L+1),       down(i)=i+L; else down(i)=i-L; end
        end
    case 'Dirichlet'
        for i=1:L^2
                if(mod(i,L)==1 || mod(i,L)==0 || i>(L-1)*L || i<L+1)
                    left(i)=i;right(i)=i;up(i)=i;down(i)=i;
                else
                    left(i)=i-1;right(i)=i+1; up(i)=i+L;down(i)=i-L;
                end
        end
end
s.left=left;
s.right=right;
s.up=up;
s.down=down;

%initial conditions for ode15s
v=vertcat(n,a,m,c,g);            
tspan=[0 tmax];
s.L=L;

s.rn=s.Dn/(dx*dx);
s.rm=s.Dm/(dx*dx);
s.rc=s.Dc/(dx*dx);
s.rg=s.Dg/(dx*dx);

% For flux limiter code
s.chem_n=s.theta_n/(dx);
s.chem_m=s.theta_m/(dx);
s.dx=dx;

options=odeset('Vectorized','off','JPattern',createJacobianMatrixChemo(L),'RelTol',1e-8,'AbsTol',1e-9);

% Initialise results matrix
w=zeros(numel(v),numTimesteps);
w(:,1)=v;

% Timestepping loop
for i=1:numTimesteps-1
    
    disp(t(i));
   
    % Explicit solve for chemotaxis terms on [t_i, t_i + dt/2]
    k1=(dt/2)*mod3_chem2d_fluxLimiter_explicit(t(i),w(:,i),s);
    k2=(dt/2)*mod3_chem2d_fluxLimiter_explicit(t(i)+dt/4,w(:,i)+k1/2,s);
    k3=(dt/2)*mod3_chem2d_fluxLimiter_explicit(t(i)+dt/4,w(:,i)+k2/2,s);
    k4=(dt/2)*mod3_chem2d_fluxLimiter_explicit(t(i)+dt/2,w(:,i)+k3,s);
    z1=w(:,i)+(1/6)*(k1+2*k2+2*k3+k4);
    
    % Implicit solve for non-chemotaxis terms on [t_i, t_i + dt]
    [~,dv]=ode15s(@(t,v)mod3_chem2d_fluxLimiter_implicit(t,v,s),[t(i) t(i+1)],z1,options);
    z2=dv(end,:)';
    
    % Explicit solve for chemotaxis terms on [t_i + dt/2, t_i + dt]
    k1=(dt/2)*mod3_chem2d_fluxLimiter_explicit(t(i)+dt/2,z2,s);
    k2=(dt/2)*mod3_chem2d_fluxLimiter_explicit(t(i)+dt/2+dt/4,z2+k1/2,s);
    k3=(dt/2)*mod3_chem2d_fluxLimiter_explicit(t(i)+dt/2+dt/4,z2+k2/2,s);
    k4=(dt/2)*mod3_chem2d_fluxLimiter_explicit(t(i)+dt/2+dt/2,z2+k3,s);
    w(:,i+1)=z2+(1/6)*(k1+2*k2+2*k3+k4);
    
    % If system has gone to zero, stop
    if(norm(w(:,i+1),Inf)<1e-8)
        break;
    end
end

w=w';
n1=w(:,1:L^2);
a1=w(:,L^2+1:2*L^2);
m1=w(:,2*L^2+1:3*L^2);
c1=w(:,3*L^2+1:4*L^2);
g1=w(:,4*L^2+1:5*L^2);

% Windows...
save(strcat('E:\Anahita\DataFiles\ParamSet2\',fname),'n1','a1','m1','c1','g1','t','Lx','-v7.3');
% Linux...
%save(strcat('/data/phy3nelsom/AnahitaDataFiles/ParamSet2/',fname),'n1','a1','m1','c1','g1','t','Lx','-v7.3');
%profile viewer
end
