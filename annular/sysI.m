function [taux] = sysI(a, m, zeta, J, Rei, Ni, kk)
% system I from "Lubricated pipelining: stability of core annular flow"
Reo = (zeta / m) * Rei;
Lo = (a - 1.0);
No = round( (a - 1.0) * Ni ); % number of inner points in outer region
n = 0; % BC are only for axisymmetric
addpath('../Chebyshev')

% we only impose the B.C. at r=a in all matrices
% the variable have 2 DOF at r=1
%%% Matrices (inner region): %%%
[ri, DMi] = chebdif(Ni+2, 3); % using an additional point at each boundary
ri = 0.5 * (ri + 1.0); % rescaling interval
for l=1:3
 DMi(:,:,l) *= (2.0/1.0)^l; % rescaling derivative
end
D0i = eye(Ni+2);
D1i=DMi(:,:,1);D2i=DMi(:,:,2);D3i=DMi(:,:,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Matrices (outer region): %%%
[ro, DMo] = chebdif(No+2, 3); % using an additional point at each boundary
DMo = DMo(2:end,2:end,:); % u(a) = v(a) = w(a) = 0
ro = ro(2:end);
ro = Lo * 0.5 * (ro + 1.0) + 1.0; % rescaling interval
for l=1:3
 DMo(:,:,l) *= (2.0/Lo)^l; % rescaling derivative
end
D0o = eye(No+1);
D1o=DMo(:,:,1);D2o=DMo(:,:,2);D3o=DMo(:,:,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = [ro;ri(2:end)]; % length(r) = (No + 1) + (Ni + 2 -1)
Oi = zeros(No+1,Ni+2);
Oo = zeros(Ni+2,No+1);
%D0 = eye(Ni+No+3);
ii = No + 1; %one = r(ii)

%%% Base flow: %%%
[Wi, dWi] = baseFlowI(ri,a,m);
[Wo, dWo] = baseFlowI(ro,a,m);
[W, dW] = baseFlowI(r,a,m);
[Wat1, dW2at1] = baseFlowI( 1.0, a, m);
dWjump = 2.0*(1.0-m)/(a-(1.0-m)/a);
d2Wi = -2.0 * m / (a - (1.0-m)/a);
d2Wo = -2.0 / (a - (1.0-m)/a);
%%%%%%%%%%%%%%%%%%

for ik = 1:length(kk)
 k = kk(ik);

 %%% Matrices of Eq. (5.12): %%%
 Mi = (1i/Rei) * ( D2i + diag(1.0./ri) * D1i - k^2 * D0i - (n^2 + 1.0) * diag(ri.^(-2)) ) + ...
  diag(dWi/k) * D1i + (d2Wi/k) * D0i + k * diag(Wi);
 Mo = (1i/Reo) * ( D2o + diag(1.0./ro) * D1o - k^2 * D0o - (n^2 + 1.0) * diag(ro.^(-2)) ) + ...
  diag(dWo/k) * D1o + (d2Wo/k) * D0o + k * diag(Wo);
 Auu = [[Mo Oi]; ...
        [Oo Mi]];

 Mi = -2i * (n/Rei) * diag(ri.^(-2));
 Mo = -2i * (n/Reo) * diag(ro.^(-2));
 Auv = [[Mo Oi]; ...
        [Oo Mi]];

 Mi = (1i/(k*Rei)) * ( D3i + diag(1.0./ri) * D2i - k^2 * D1i - (n^2 + 1.0) * diag(ri.^(-2)) * D1i + 2.0 * n^2 * diag(ri.^(-3)) ) + ...
  diag(Wi) * D1i + diag(dWi);
 Mo = (1i/(k*Reo)) * ( D3o + diag(1.0./ro) * D2o - k^2 * D1o - (n^2 + 1.0) * diag(ro.^(-2)) * D1o + 2.0 * n^2 * diag(ro.^(-3)) ) + ...
  diag(Wo) * D1o + diag(dWo);
 Auw = [[Mo Oi]; ...
        [Oo Mi]];

 Mi = k * D0i;
 Mo = k * D0o;
 Buu = [[Mo Oi]; ...
        [Oo Mi]];

 Buv = zeros(size(Buu));

 Mi = D1i;
 Mo = D1o;
 Buw = [[Mo Oi]; ...
        [Oo Mi]];

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% Matrices of Eq. (5.13): %%%
 Mi = ( (2i/Rei) * diag(ri.^(-2)) + diag(dWi/k) * diag(1.0./ri) ) * n;
 Mo = ( (2i/Reo) * diag(ro.^(-2)) + diag(dWo/k) * diag(1.0./ro) ) * n;
 Avu = [[Mo Oi]; ...
        [Oo Mi]];

 Mi = (-1i/Rei) * ( D2i + diag(1.0./ri) * D1i - k^2 * D0i - (n^2 + 1.0) * diag(ri.^(-2)) ) - k * diag(Wi);
 Mo = (-1i/Reo) * ( D2o + diag(1.0./ro) * D1o - k^2 * D0o - (n^2 + 1.0) * diag(ro.^(-2)) ) - k * diag(Wo);
 Avv = [[Mo Oi]; ...
        [Oo Mi]];

 Mi = (n*1i/(k*Rei)) * ( diag(1.0./ri) * D2i + diag(ri.^(-2)) * D1i - k^2 * diag(1.0./ri) - n^2 * diag(ri.^(-3)) ) + diag(Wi);
 Mo = (n*1i/(k*Reo)) * ( diag(1.0./ro) * D2o + diag(ro.^(-2)) * D1o - k^2 * diag(1.0./ro) - n^2 * diag(ro.^(-3)) ) + diag(Wo);
 Avw = [[Mo Oi]; ...
        [Oo Mi]];

 Bvu = zeros(size(Buu));

 Mi = -k * D0i;
 Mo = -k * D0o;
 Bvv = [[Mo Oi]; ...
        [Oo Mi]];

 Mi = n * diag(1.0./ri);
 Mo = n * diag(1.0./ro);
 Bvw = [[Mo Oi]; ...
        [Oo Mi]];

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% Matrices of Eq. (5.2): %%%
 Mi = D1i + diag(1.0./ri);
 Mo = D1o + diag(1.0./ro);
 Awu = [[Mo Oi]; ...
        [Oo Mi]];

 Mi = diag(n./ri);
 Mo = diag(n./ro);
 Awv = [[Mo Oi]; ...
        [Oo Mi]];

 Mi = k * D0i;
 Mo = k * D0o;
 Aww = [[Mo Oi]; ...
        [Oo Mi]];

 Bwu = zeros(size(Awu));
 Bwv = zeros(size(Awu));
 Bww = zeros(size(Awu));
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 A = [ [Auu Auv Auw]; ...
       [Avu Avv Avw]; ...
       [Awu Awv Aww] ];

 B = [ [Buu Buv Buw]; ...
       [Bvu Bvv Bvw]; ...
       [Bwu Bwv Bww] ];

 %%% solving A u = B c u: %%%
 [eve,eva] = eig(A,B);
 c = diag(eva);
 omega = c * k;
 tau = imag(omega);
 tau = tau .* isfinite(tau);
 tau = tau .* (tau < 0.03);
 [taux(ik), ip] = max(tau);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% plotting eigenvectors: %%%
% lgth = size(eve,1) / 3;
% u = eve(1:lgth,ip);
% u = u([1:ii, ii+2:lgth]);
% figure(1)
% plot(r, u,'bk'); hold on; plot(r, W, 'r')
% ylabel('w'); xlabel('r');legend('perturbation','base flow')
% pause(1)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

