function p=setParameterValues_ParamSet2(varargin)

%pro-inflammatory mediators (c) parameters
p.alpha=0;                     %production rate [/]
p.nu=0.1;                      % apoptosis rate [/]

%active neutrophils (n) parameters
p.kn=0.01;                      %neutrophils mediators production rate [/]
p.beta_n=0.1;                    %saturation constant [/]

%apoptotic neutrophils (a) parameters
p.gamma_a=1;                    %necrosis rate [/]
p.phi=0.1;                       %removal rate [/]
p.beta_a=0.1;                    %saturation constant [/]  

%macrophages (m) parameters
p.gamma_m=0.01;               %leaving tissue rate [/]

%anti-inflammatory mediators (g) parameters
p.beta_c=1.2e-1;                      %saturation constant [/]
p.beta_g=0.01;
p.gamma_g=1;                   %decay rate [/]
p.kg=0.1;                      %production rate [/]

% Spatial parameters
p.Dc=1e-4;
p.Dg=1e-4;
p.Dn=1e-5;
p.Dm=1e-6;
p.theta_n=1e-5; 
p.theta_m=1e-6;

if(nargin>0)
   if(strcmp(varargin{1},'phi'))
       p.phi=varargin{2};
   end
   if(strcmp(varargin{1},'nu'))
       p.nu=varargin{2};
   end 
   if(nargin>2)
       if(strcmp(varargin{3},'phi'))
           p.phi=varargin{4};
       end
       if(strcmp(varargin{3},'nu'))
           p.nu=varargin{4};
       end    
   end
end


% y0=[0.310848,0,21.194060,0.262276,0.0730]';
% y=fsolve(@(y)model3_odes(0,y,p,0),y0);
% n0=y(1);
% a0=y(2);
% m0=y(3);
% c0=y(4);
% g0=y(5);

end