clear
close all

%% Parameter settings
% Domain parameters
nx=308; % domain size along x
ny=218; % domain size along y

% Well parameters
x_well=nx/2; % pumping well location along x
y_well=ny/2; % pumping well location along y
Q=4*pi; % pumping rate
T=1; % transmissivity (assumed unity)
S=5e-4; % storativity
rw=0.05; % well screen radius (for finite well setting)
rc=0.05; % well casing radius (for finite well setting)

% Parameters for gradient interpolation
lam_topo=30; % correlation scale for gradient interpolation
err_topo=1e-2; % error variance for gradient interpolation

% Parameters for topography interpolation
lam_est=40; % correlation scale for topographic interpolation
err_est=1e-1; % error variance for topographic interpolation

%% Loading point data (aquifer heights & orientations)
load_data

%% Gradient interpolation
% Interpolation using grad_dat (pointwise aquifer orientation)
dfdx_grid=GPR_est(0,nx,ny,grad_pnt(:,[1 2 3]),lam_topo,err_topo,[]);
dfdy_grid=GPR_est(0,nx,ny,grad_pnt(:,[1 2 4]),lam_topo,err_topo,[]);
df_grid=cat(3,dfdx_grid, dfdy_grid);

%% Aquifer topography interpolation
% Interpolation using dat_trn and gradient information
% Running geodesic-kernel-based GPR (Piao and Park, 2023)
T_est=GPR_est(1,nx,ny,dat_pnt,lam_est,err_est,df_grid);
T_est=reshape(T_est,ny,nx); % reshape topgraphic vector into nxXny
[dTdx,dTdy]=gradient(T_est); % compute topographic gradients

%% Computing geodesic distances from pumping well
[xx,yy]=meshgrid(1:nx,1:ny);
xx=xx(:);
yy=yy(:);

dg_w=max(0.1,comp_d_g_aniso(nx,ny,x_well,y_well,xx,yy,dTdx,dTdy,1,1,0));
% For anisotropy, you may use this line instead of the above line of dg_w
% dg_w=max(0.1,comp_d_g_aniso(nx,ny,x_well,y_well,xx,yy,dTdx,dTdy,2,1,pi/4));
% where scale factor along x and y are 2 and 1, respectively, and the
% rotation angle is pi/4A
%% Computing drawdowns
ddn=zeros(nx*ny,1);
for ii=1:nx*ny
    f=basic_wh(2);
    ddn(ii)=talbot_Lap_inv(@(p) f(dg_w(ii),p,Q,S,rw,rc),1);
end
ddn=reshape(ddn,ny,nx);


%% Drawing Figure for estimated Z
draw_figures

%% Well-hydraulics solutions
% References
%   Agarwal et al. (1970)
%   Hantush (1960)
%   Hantush and Jacob (1955)
%   Papadoupulos and Cooper (1967)
%   Theis (1935)
% cs: analytical solution identifier
% f_wh: drawdowns based on the selected analytical solution
function f_wh=basic_wh(cs)
switch cs
    case 1 % Theis
        f_wh=@(d_g,p,Q,S) Q/(2*pi*p)*besselk(0,d_g*sqrt(S*p));
    case 2 % finite diameter (wellbore storage)
        f_wh=@(d_g,p,Q,S,rw,rc) Q/(2*pi*p)*besselk(0,d_g*sqrt(S*p))/...
            (rw*sqrt(S*p)*besselk(1,rw*sqrt(S*p))+...
            rc^2/(2)*p*besselk(0,rw*sqrt(S*p)));
    case 3 % a finite diameter well with skin
        f_wh=@(d_g,p,Q,S,rw,rc) Q/(2*pi*p)*(besselk(0,d_g*sqrt(S*p))+...
            sig*rw*sqrt(S*p)*besselk(1,rw*sqrt(S*p)));
    case 4 % a leaky confined aquifer with infinitesimal storage
        f_wh=@(d_g,p,Q,S,B) Q/(2*pi*p)*...
            besselk(0,d_g*sqrt(S*p+1/B^2));
    case 5 % a leaky confined aquifer with finite storage
        f_wh=@(d_g,p,Q,S,Kl,Ssl,bl) Q/(2*pi*p)*...
            besselk(0,d_g*sqrt(S*p+Kl*sqrt(Ssl/Kl*p)*coth(sqrt(Ssl/Kl*p)*bl)));
    case 6 % a dual porosity confined aquifer
        f_wh=@(d_g,p,Q,Ss,b,Ssm,alp) Q/(2*pi*b*p)*...
            besselk(0,d_g*sqrt(Ss*p+alp*Ssm*p/(alp+Ssm*p)));
end
end

%% Talbot numerical inverse Laplace transformation
% f_s: Laplace transformed function handle
% t: time points vector where the inverse transform is evaluated
% f_t: inverse Laplace transform values at time points t
function f_t=talbot_Lap_inv(f_s,t)
M=16; % Number of terms in the series
k=1:(M-1);  % Iteration index
% Calculate delta for every index
delta=zeros(1,M);
delta(1)=2 * M / 5;
delta(2:end)=2*pi/5*k.*(cot(pi/M*k)+1i);
% Calculate gamma for every index
gamma=zeros(1,M);
gamma(1)=0.5*exp(delta(1));
gamma(2:end)=(1+1i*pi/M*k.*(1+cot(pi/M*k).^2)-...
             1i*cot(pi/M*k)).*exp(delta(2:end));
% Create meshgrids for delta and t
[delta_mesh,t_mesh]=meshgrid(delta,t);
gamma_mesh=meshgrid(gamma,t);
% Calculate inverse Laplace transform
f_t=0.4./t.*sum(real(gamma_mesh.*...
                     arrayfun(@(z) f_s(z),delta_mesh./t_mesh)),2);
end
