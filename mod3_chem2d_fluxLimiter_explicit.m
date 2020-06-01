%model 3 space 2D discretized system equation (ODEs in time) 

function dv=mod3_chem2d_fluxLimiter_explicit(t,v,s)

n=v(1:s.L^2);
%a=v(s.L^2+1:2*s.L^2);
m=v(2*s.L^2+1:3*s.L^2);
c=v(3*s.L^2+1:4*s.L^2);
%g=v(4*s.L^2+1:5*s.L^2);

phi=@(r)(abs(r)+r)./(1+abs(r));

nfluxH=(c(s.right)-c(s.left)).*n/(2*s.dx);
nfluxV=(c(s.up)-c(s.down)).*n/(2*s.dx);
mfluxH=(c(s.right)-c(s.left)).*m/(2*s.dx);
mfluxV=(c(s.up)-c(s.down)).*m/(2*s.dx);

eps=1e-16;
nrH=(nfluxH(s.right)-nfluxH + eps)./(nfluxH-nfluxH(s.left) + eps);
nrV=(nfluxV(s.up)-nfluxV + eps)./(nfluxV-nfluxV(s.down) + eps);
mrH=(mfluxH(s.right)-mfluxH + eps)./(mfluxH-mfluxH(s.left) + eps);
mrV=(mfluxV(s.up)-mfluxV + eps)./(mfluxV-mfluxV(s.down) + eps);

nSwitchH_pos=((n>0 & nfluxH>0) | (n==0 & nfluxH(s.left)>0));
nSwitchH_neg=((n>0 & nfluxH<0) | (n==0 & nfluxH(s.right)<=0));
nSwitchV_pos=((n>0 & nfluxV>0) | (n==0 & nfluxV(s.down)>0));
nSwitchV_neg=((n>0 & nfluxV<0) | (n==0 & nfluxV(s.up)<=0));

mSwitchH_pos=((m>0 & mfluxH>0) | (m==0 & mfluxH(s.left)>0));
mSwitchH_neg=((m>0 & mfluxH<0) | (m==0 & mfluxH(s.right)<=0));
mSwitchV_pos=((m>0 & mfluxV>0) | (m==0 & mfluxV(s.down)>0));
mSwitchV_neg=((m>0 & mfluxV<0) | (m==0 & mfluxV(s.up)<=0));

nchemc=(nSwitchH_pos).*(1+0.5*phi(nrH)-phi(nrH(s.left))./(2*nrH(s.left))).*(nfluxH-nfluxH(s.left))+...
       (nSwitchH_neg).*(1+0.5*phi(1./nrH)-phi(1./nrH(s.right))./(2./nrH(s.right))).*(nfluxH(s.right)-nfluxH)+...
       (nSwitchV_pos).*(1+0.5*phi(nrV)-phi(nrV(s.down))./(2*nrV(s.down))).*(nfluxV-nfluxV(s.down))+...
       (nSwitchV_neg).*(1+0.5*phi(1./nrV)-phi(1./nrV(s.up))./(2./nrV(s.up))).*(nfluxV(s.up)-nfluxV);
   
mchemc=(mSwitchH_pos).*(1+0.5*phi(mrH)-phi(mrH(s.left))./(2*mrH(s.left))).*(mfluxH-mfluxH(s.left))+...
       (mSwitchH_neg).*(1+0.5*phi(1./mrH)-phi(1./mrH(s.right))./(2./mrH(s.right))).*(mfluxH(s.right)-mfluxH)+...
       (mSwitchV_pos).*(1+0.5*phi(mrV)-phi(mrV(s.down))./(2*mrV(s.down))).*(mfluxV-mfluxV(s.down))+...
       (mSwitchV_neg).*(1+0.5*phi(1./mrV)-phi(1./mrV(s.up))./(2./mrV(s.up))).*(mfluxV(s.up)-mfluxV);   

dv(1:s.L^2)=-s.chem_n*nchemc;

dv(s.L^2+1:2*s.L^2)=0;

dv(2*s.L^2+1:3*s.L^2)=-s.chem_m*mchemc;

dv(3*s.L^2+1:4*s.L^2)=0;

dv(4*s.L^2+1:5*s.L^2)=0;

dv=dv';
