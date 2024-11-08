% Computing geodesic distance
% nx, ny: domain size
% u0, v0: current location
% u1, v1: locations to be computed for geodesic distances from (u0, v0)
% df_du, df_dv: manifold gradient information
% alp, bet: anisotropy factors for the u- and v-directions, respectively
% theta: rotation angle of the anisotropy
function d_g = comp_d_g_aniso(nx, ny, u0, v0, u1, v1, df_du, df_dv, alp, bet, theta)

nd = 50; % number of segments

% Legendre-Gauss Quadrature Weights and Nodes
[absc, wght] = lgwt(nd, -1, 1);
ulam = @(t) u0 + (u1 - u0) * t'; % Lambda function for u(lambda)
vlam = @(t) v0 + (v1 - v0) * t'; % Lambda function for v(lambda)
t = 0.5 * (absc + 1); % Scaling for quadrature points

ult = ulam(t); % u(lambda)
vlt = vlam(t); % v(lambda)

% Interpolating manifold gradient information at u(lambda), v(lambda)
dfdu = reshape(df_du, ny, nx);
dfdv = reshape(df_dv, ny, nx);
df_du_loc = interp2(reshape(u1, ny, nx), reshape(v1, ny, nx), dfdu, ult(:), vlt(:));
df_dv_loc = interp2(reshape(u1, ny, nx), reshape(v1, ny, nx), dfdv, ult(:), vlt(:));

% Adjusting for gradient normalization (if needed)
df_du_loc = -1 * df_du_loc;
df_dv_loc = -1 * df_dv_loc;

% Compute components of the geometric metric tensor
g11 = 1 + reshape(df_du_loc, nx * ny, nd).^2;
g12 = reshape(df_du_loc, nx * ny, nd) .* reshape(df_dv_loc, nx * ny, nd);
g22 = 1 + reshape(df_dv_loc, nx * ny, nd).^2;

% Anisotropy matrix
A = [1 / alp, 0; 0, 1 / bet];

% Rotation matrix with an angle theta
R = [cos(theta), -sin(theta); sin(theta), cos(theta)];

% Apply rotation to the anisotropy matrix
A_rot = R * A * R'; % Rotated anisotropy matrix

% Pre-allocate new matrices for g11, g12, g22 after transformation
g11_transformed = zeros(size(g11));
g12_transformed = zeros(size(g12));
g22_transformed = zeros(size(g22));

% Use parfor for parallel computation of element-wise transformations
parfor i = 1:(nx * ny)
    for j = 1:nd
        % Build the geometric metric tensor at this grid point
        g_local = [g11(i, j), g12(i, j); g12(i, j), g22(i, j)];
        
        % Apply anisotropy and rotation transformations
        g_local_transformed = A_rot * g_local * A_rot';
        
        % Extract transformed values
        g11_transformed(i, j) = g_local_transformed(1, 1);
        g12_transformed(i, j) = g_local_transformed(1, 2);
        g22_transformed(i, j) = g_local_transformed(2, 2);
    end
end

% Integration of geodesic distance using Legendre-Gauss quadrature
d_g = sqrt((ult - u0).^2 .* g11_transformed + (vlt - v0).^2 .* g22_transformed + ...
            2 * (ult - u0) .* (vlt - v0) .* g12_transformed) * 0.5 * wght;

% Sum over all segments to get the final geodesic distance
d_g = sum(d_g, 2);

end

%% Legendre-Gauss abscissa and weights
% N: number of segments
% a, b: range
% x: abscissa
% w: weights
function [x,w]=lgwt(N,a,b)

N=N-1;
N1=N+1; N2=N+2;
xu=linspace(-1,1,N1)';
% Initial guess
y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);
% Legendre-Gauss Vandermonde Matrix
L=zeros(N1,N2);
% Derivative of LGVM
Lp=zeros(N1,N2);
% Compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method
y0=2;
% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0))>eps
    L(:,1)=1;
    Lp(:,1)=0;
    L(:,2)=y;
    Lp(:,2)=1;
    for k=2:N1
        L(:,k+1)=((2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1))/k;
    end
    Lp=(N2)*(L(:,N1)-y.*L(:,N2))./(1-y.^2);   
    y0=y;
    y=y0-L(:,N2)./Lp;
end
% Linear map from[-1,1] to [a,b]
x=(a*(1-y)+b*(1+y))/2;      
% Compute the weights
w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
end
