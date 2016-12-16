function [taux, cx] = sysI(n, a, m, zeta, J, Rei, Ni, kk)
% system I from "Lubricated pipelining: stability of core annular flow"
Reo = (zeta / m) * Rei;
Lo = (a - 1.0);
No = round( (a - 1.0) * Ni ); % number of inner points in outer region
addpath('../Chebyshev')

% we only impose the B.C. at r=a in all matrices
% the variable have 2 DOF at r=1
%%% Matrices (inner region): %%%
[ri, DMi] = chebdif(Ni+2, 3); % using an additional point at each boundary
nn = 0;
if ( n > 1 )
 nn = 1;
 DMi = DMi(1:end-1,1:end-1,:); % u(0) = v(0) = w(0) = 0
 ri = ri(1:end-1);
end
ri = 0.5 * (ri + 1.0); % rescaling interval
for l=1:3
 DMi(:,:,l) *= (2.0/1.0)^l; % rescaling derivative
end
D0i = eye(Ni+2-nn);
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
Oi = zeros(No+1,Ni+2-nn);
Oo = zeros(Ni+2-nn,No+1);
ii = No + 1; %one = r(ii)

%%% Base flow: %%%
[Wi, dWi] = baseFlowI(ri,a,m);
[Wo, dWo] = baseFlowI(ro,a,m);
[W, dW] = baseFlowI(r,a,m);
[Wat1, dW2at1] = baseFlowI( 1.0, a, m);
dWjump = 2.0*(1.0-m)/(a*a-(1.0-m));
dWzetaJump = 2.0*(1.0-m*zeta)/(a*a-(1.0-m));
d2Wi = -2.0 * m / (a*a - (1.0-m));
d2Wo = -2.0 / (a*a - (1.0-m));
%%%%%%%%%%%%%%%%%%

for ik = 1:length(kk)
 k = kk(ik);

 %%% Matrices of Eq. (5.12): %%%
 Mi = (1i/Rei) * ( D2i + diag(1.0./ri) * D1i - k^2 * D0i - (n^2 + 1.0) * diag(ri.^(-2)) ) + diag(dWi/k) * D1i + (d2Wi/k) * D0i + k * diag(Wi);
 Mo = (1i/Reo) * ( D2o + diag(1.0./ro) * D1o - k^2 * D0o - (n^2 + 1.0) * diag(ro.^(-2)) ) + diag(dWo/k) * D1o + (d2Wo/k) * D0o + k * diag(Wo);
 Auu = [[Mo Oi]; ...
        [Oo Mi]];

 Mi = -2i * (n/Rei) * diag(ri.^(-2));
 Mo = -2i * (n/Reo) * diag(ro.^(-2));
 Auv = [[Mo Oi]; ...
        [Oo Mi]];

 Mi = (1i/(k*Rei)) * ( D3i + diag(1.0./ri) * D2i - k^2 * D1i - (n^2 + 1.0) * diag(ri.^(-2)) * D1i + 2.0 * n^2 * diag(ri.^(-3)) ) + diag(Wi) * D1i + diag(dWi);
 Mo = (1i/(k*Reo)) * ( D3o + diag(1.0./ro) * D2o - k^2 * D1o - (n^2 + 1.0) * diag(ro.^(-2)) * D1o + 2.0 * n^2 * diag(ro.^(-3)) ) + diag(Wo) * D1o + diag(dWo);
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

 Mi = (n*1i/(k*Rei)) * ( diag(1.0./ri) * D2i + diag(ri.^(-2)) * D1i - k^2 * diag(1.0./ri) - n^2 * diag(ri.^(-3)) ) + n * diag(1.0./ri) * diag(Wi);
 Mo = (n*1i/(k*Reo)) * ( diag(1.0./ro) * D2o + diag(ro.^(-2)) * D1o - k^2 * diag(1.0./ro) - n^2 * diag(ro.^(-3)) ) + n * diag(1.0./ro) * diag(Wo);
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
 %%% B.C. at r = 0 %%%
 if(n==0)
 Auu(end,:) = 0;Auu(end,end) = 1;
 Auv(end,:) = 0;
 Auw(end,:) = 0;
 Avu(end,:) = 0;
 Avv(end,:) = 0;Avv(end,end) = 1;
 Avw(end,:) = 0;
 Awu(end,:) = 0;
 Awv(end,:) = 0;
 Aww(end,:) = 0;Aww(end,end) = k;
 Awu(end,end-size(D1i,2)+1:end) = 2.0 * D1i(end,:);

 Buu(end,:) = 0;
 Buv(end,:) = 0;
 Buw(end,:) = 0;
 Bvu(end,:) = 0;
 Bvv(end,:) = 0;
 Bvw(end,:) = 0;
 Bwu(end,:) = 0;
 Bwv(end,:) = 0;
 Bww(end,:) = 0;
 elseif (n == 1)
 return
 end
 %%%%%%%%%%%%%%%%%%%%%
 %%% B.C. (5.7-8) %%%
 Auu(ii,:) = 0;
 Auv(ii,:) = 0;
 Auw(ii,:) = 0;
 Avu(ii,:) = 0;
 Avv(ii,:) = 0;
 Avw(ii,:) = 0;
 Awu(ii,:) = 0;
 Awv(ii,:) = 0;
 Aww(ii,:) = 0;
 Buu(ii,:) = 0;
 Buv(ii,:) = 0;
 Buw(ii,:) = 0;
 Bvu(ii,:) = 0;
 Bvv(ii,:) = 0;
 Bvw(ii,:) = 0;
 Bwu(ii,:) = 0;
 Bwv(ii,:) = 0;
 Bww(ii,:) = 0;
 Auu(ii,ii) = 1;Auu(ii,ii+1) = -1;
 Avv(ii,ii) = 1;Avv(ii,ii+1) = -1;
 Awu(ii,ii) = dWjump;
 Aww(ii,ii) = -Wat1 * k;Aww(ii,ii+1) = Wat1 * k;
 Bww(ii,ii) = -k;Bww(ii,ii+1) = k;
 %%%%%%%%%%%%%%%%%%%%
 %%% B.C. (5.9-10) & (5.14) %%%
 Auu(ii+1,:) = 0;
 Auv(ii+1,:) = 0;
 Auw(ii+1,:) = 0;
 Avu(ii+1,:) = 0;
 Avv(ii+1,:) = 0;
 Avw(ii+1,:) = 0;
 Awu(ii+1,:) = 0;
 Awv(ii+1,:) = 0;
 Aww(ii+1,:) = 0;
 Buu(ii+1,:) = 0;
 Buv(ii+1,:) = 0;
 Buw(ii+1,:) = 0;
 Bvu(ii+1,:) = 0;
 Bvv(ii+1,:) = 0;
 Bvw(ii+1,:) = 0;
 Bwu(ii+1,:) = 0;
 Bwv(ii+1,:) = 0;
 Bww(ii+1,:) = 0;
 %% (5.9) %%
 Auw(ii+1,1:ii) = -m * D1o(ii,:);
 Auw(ii+1,ii+1:end) = D1i(1,:);
 Auu(ii+1,1:ii) = m * k * D0o(ii,:);
 Auu(ii+1,ii+1:end) = -k * D0i(1,:);
 %%%%%%%%%%%
 %% (5.10) %%
 Avv(ii+1,1:ii) = -m * (D1o(ii,:) - D0o(ii,:));
 Avv(ii+1,ii+1:end) = D1i(1,:) - D0i(1,:);
 Avu(ii+1,1:ii) = m * n * D0o(ii,:);
 Avu(ii+1,ii+1:end) = -n * D0i(1,:);
 %%%%%%%%%%%%
 Aww(ii+1,1:ii) = -( J *(1.0 - k^2 - n^2) * D0o(ii,:) / (dWjump * Rei^2) + m * (1i/k) * (D2o(ii,:) + D1o(ii,:) - (k^2 + n^2) * D0o(ii,:)) + Wat1 * (zeta - dWzetaJump / dWjump) * D0o(ii,:) );
 Aww(ii+1,ii+1:end) = J *(1.0 - k^2 - n^2) * D0i(1,:) / (dWjump * Rei^2) + (1i/k) * (D2i(1,:) + D1i(1,:) - (k^2 + n^2) * D0i(1,:)) + Wat1 * (1.0 - dWzetaJump / dWjump) * D0i(1,:);
 Awu(ii+1,1:ii) = -2i * m  * D1o(ii,:);
 Awu(ii+1,ii+1:end) = 2i * D1i(1,:);
 Bww(ii+1,1:ii) =   -(zeta - dWzetaJump / dWjump) * D0o(ii,:);
 Bww(ii+1,ii+1:end) = (1.0 - dWzetaJump / dWjump) * D0i(1,:);
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
 tau = tau .* (tau < 1.3);
 [taux(ik), ip] = max(tau);
 cx = c(ip);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% plotting eigenvectors: %%%
 %lgth = size(eve,1) / 3;
 %u = eve(1:lgth,ip);
 %u(ii) - u(ii+1)
% u = u([1:ii, ii+2:lgth]);
% figure(1)
% plot(r, u,'bk'); hold on; plot(r, W, 'r')
% ylabel('w'); xlabel('r');legend('perturbation','base flow')
% pause(1)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

