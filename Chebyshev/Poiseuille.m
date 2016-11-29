% Solving Poiseuille flow with Chebyshev
px = -1.4; % pressure gradient
a = 1.3; % outer radius (inner radius = 1)
m = 2.2; % viscosity ratio

Ni = 600; % number interior points
No = 50; % number exterior points

mu = 1.0/m;
Lo = (a - 1.0);
Oo = zeros(Ni,No);
Oi = zeros(No,Ni);
[ri, DM] = chebdif(Ni+2, 4);
DM(:,end-1,:) += DM(:,end,:); % Neumann BC
DM = DM(1:end-1,1:end-1,:);
ri = ri(1:end-1);
% rescaling interval
ri = 0.5 * (ri + 1.0);
for n=1:4
 DM(:,:,n) *= (2.0/1.0)^n;
end
D1i=DM(:,:,1);D2i=DM(:,:,2);D3i=DM(:,:,3);D4i=DM(:,:,4);
[ro, DM] = chebdif(No+2, 4);
DM = DM(2:end,2:end,:); % no slip BC
ro = ro(2:end-1);
% rescaling interval
ro = Lo * 0.5 * (ro + 1.0) + 1.0;
for n=1:4 % rescaling:
 DM(:,:,n) *= (2.0/Lo)^n;
end
D1o=DM(:,:,1);D2o=DM(:,:,2);D3o=DM(:,:,3);D4o=DM(:,:,4);

Ai = mu * (diag(ri(2:end)) * D2i(2:end,:) + D1i(2:end,:));
Ao =      (diag(ro) * D2o(1:end-1,:) + D1o(1:end-1,:));
% imposing interface BC with iBC:
iBC = zeros(1, Ni+No+1);
iBC(1:No+1) += D1o(end,:);
iBC(No+1:end) -= mu * D1i(1,:);
A = [[Ao Oi];iBC;[Oo Ai]];
r = [ro;ri];
%b = px * ones(length(r),1);
b = px * r;
b(No+1) = 0;
u = A \ b;
addpath('../annular')
ua = baseFlowI(r,a,m);
plot(r, u, r, ua,'g')
