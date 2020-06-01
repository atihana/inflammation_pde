function dydt=model3_odes(t,y,p,f)

n=y(1);
a=y(2);
m=y(3);
c=y(4);
g=y(5);

dydt(1)=c./(1+g) - p.nu*n.*(1+g/p.beta_g)./(1+c/p.beta_c);
dydt(2)=p.nu*n.*(1+g/p.beta_g)./(1+c/p.beta_c)-p.gamma_a*a-p.phi*m.*a;
dydt(3)=c-p.gamma_m*m;
%dydt(4)=p.alpha*f(t,p.A)+p.kn*(n.^2./(p.beta_n^2+n.^2))+p.gamma_a*(a.^2./(p.beta_a^2+a.^2))-c;
dydt(4)=p.kn*(n.^2./(p.beta_n^2+n.^2))+p.gamma_a*(a.^2./(p.beta_a^2+a.^2))-c;
dydt(5)=p.kg*p.phi*m.*a-p.gamma_g*g;

dydt=dydt';

end
