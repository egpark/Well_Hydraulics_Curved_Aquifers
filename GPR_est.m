%% Talbot numerical inverse Laplace transformation
% type: identifier for conventional or new GPR (Piao and Park, 2023)
% nx,ny: domain size
% dat: pointwise data
% lam: correlation scale
% err: estimation error
% gradat: gradient data guiding estimation (Piao and Park, 2023)
% z_est: estimated values
% z_var: estimation variances
function [z_est,z_var]=GPR_est(type,nx,ny,dat,lam,err,gradat)
xo=dat(:,1); % x-coordinate
yo=dat(:,2); % y-coordinate
val=dat(:,3); % value
ndat=size(dat,1); % number of data points

idx0=(xo-1)*ny+yo;
% Computations
% Kernel matrix computations (kernel trick)
Sig_ab=ones(nx*ny,ndat+1);% Sigma_ab
Sig_bb=zeros(ndat+1);% Sigma_bb; note that the dimension is ndat '+1'
mm=0;
for ii=1:ndat
    mm=mm+1;
    Sig_ab(:,mm)=Sigma(type,idx0(ii),nx,ny,lam,gradat)';
end

Sig_bb(1:end-1,:)=Sig_ab(idx0,:);
Sig_bb(end,1:end-1)=1;
iSig_bb=pinv(Sig_bb,err);
z_est=Sig_ab*(iSig_bb*[val;0]);% value
% Relative uncertainty assessment
z_var=max(0,1-sum(Sig_ab*iSig_bb.*Sig_ab,2));
end

%% Covariance matrix computation
% type: identifier for conventional or new GPR (Piao and Park, 2023)
% idx: index for current location
% nx,ny: domain size
% lam: correlation scale
% gradat: gradient data guiding estimation (Piao and Park, 2023)
% Sig: covariance vector of idx-th column (or row)
function Sig=Sigma(type,idx,nx,ny,lam,gradat)

y0=double(mod(idx,ny));
x0=double(ceil(idx/ny));
if y0==0,y0=ny;end
[xx,yy]=meshgrid(1:nx,1:ny);
xx=xx(:);
yy=yy(:);

if type==0
    dist=sqrt((xx-x0).^2+(yy-y0).^2);
else
    dist=comp_d_g_aniso(nx,ny,x0,y0,xx,yy,gradat(:,:,1),gradat(:,:,2),1,1,0);
end

Sig=exp(-1/2*dist.^2/lam^2);
end
