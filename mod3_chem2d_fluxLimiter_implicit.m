%model 3 space 2D discretized system equation (ODEs in time) 

function dv=mod3_chem2d_fluxLimiter_implicit(t,v,s)

n=v(1:s.L^2);
a=v(s.L^2+1:2*s.L^2);
m=v(2*s.L^2+1:3*s.L^2);
c=v(3*s.L^2+1:4*s.L^2);
g=v(4*s.L^2+1:5*s.L^2);
dv(1:s.L^2)=-s.ni*((1+g/s.beta_g)./(1+c/s.beta_c)).*n+...
    c./(1+g)+s.rn*s.dmat*n;

dv(s.L^2+1:2*s.L^2)=s.ni*((1+g/s.beta_g)./(1+c/s.beta_c)).*n-...
    s.gamma_a*a-s.phi*m.*a;

dv(2*s.L^2+1:3*s.L^2)=-s.gamma_m*m+...
    c+s.rm*s.dmat*m;

dv(3*s.L^2+1:4*s.L^2)=-c+...
    s.rc*s.dmat*c+...
    s.alpha*s.dam+s.gamma_a*a.^2./(a.^2+s.beta_a^2)+...
    s.kn*n.^2./(n.^2+s.beta_n^2);

dv(4*s.L^2+1:5*s.L^2)=-s.gamma_g*g+...
    s.kg*s.phi*m.*a+s.rg*s.dmat*g;

dv=dv';
