% axisymetric system II from "Lubricated pipelining: stability of core annular flow"
Fall = 1

if (Fall == 1)
 a = 1.43;
 m = 0.5;
 J = 0;
 ReiRei = [26.42];
 kk = [0.01:0.5:25];
elseif (Fall == 2)
 a = 1.0 / 0.8; % outer radius (inner radius = 1)
 m = 0.1; % viscosity ratio
 J = 0.8*1000.0; % surface tension parameter
 ReiRei = 0.8 * [10:20:500];
 kk = [0.01:0.01:5];
end
MAXeva = 0.3; % threshold above which eigenvalues are not considered

Ni = 40; % number of points inner region
No = round( (a - 1.0) * Ni ); % number of points outer region

Fig2sig = [];
Fig2k = [];
GR = [];
for iRR = 1:length(ReiRei)
Rei = ReiRei(iRR);
Reo = (1.0/m) * Rei;
Re = Rei / 0.8; % for paper II data
Lo = (a - 1.0);

addpath('../Chebyshev')
%%% Matrices (inner region): %%%
Ii = eye(Ni);
D0i = [zeros(Ni,1) Ii];
Oi = zeros(No,Ni);
[ri, DMi] = chebdif(Ni+2, 4); % using 2 ghost points (one at each boundary)
DMi = DMi(1:end-1,1:end-1,:); % u1(0) = 0
ri = ri(2:end-1);
ri = 0.5 * (ri + 1.0); % rescaling interval
for n=1:4
 DMi(:,:,n) *= (2.0/1.0)^n; % rescaling derivative
end
D1i=DMi(:,:,1);D2i=DMi(:,:,2);D3i=DMi(:,:,3);D4i=DMi(:,:,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Matrices (outer region): %%%
Io = eye(No);
D0o = [Io zeros(No,1)];
Oo = zeros(Ni,No);
[ro, DMo] = chebdif(No+3, 4); % using 3 ghost points
% u2(a) = u2'(a)= 0:
for n=1:4
 for o=3:No+3
  for p=3:No+3
   DMo(p,o,n) -= DMo(2,o,1) .* DMo(p,2,n) / DMo(2,2,1);
  end
 end
end
DMo = DMo(3:end,3:end,:);
ro = ro(3:end-1);
ro = Lo * 0.5 * (ro + 1.0) + 1.0; % rescaling interval
for n=1:4
 DMo(:,:,n) *= (2.0/Lo)^n; % rescaling derivative
end
D1o=DMo(:,:,1);D2o=DMo(:,:,2);D3o=DMo(:,:,3);D4o=DMo(:,:,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = [ro;1.0;ri];

%%% Base flow: %%%
[Wi, dWi] = baseFlowI(ri,a,m);
[Wo, dWo] = baseFlowI(ro,a,m);
[W, dW] = baseFlowI(r,a,m);
[Wat1, dW2at1] = baseFlowI( 1.0, a, m);
%%%%%%%%%%%%%%%%%%

for ik = 1:length(kk)
 k = kk(ik);
 %%% A, B given by (5.25): %%%
 Ai = diag(ri.^4) * D4i(2:end,:) ...
      + 2.0 * diag(ri.^3) * D3i(2:end,:) ...
      - ( ( Rei * 1i * k * diag(Wi) + 2.0 * k.^2 * Ii ) * diag(ri.^2) + 3.0 * Ii ) * diag(ri.^2) * D2i(2:end,:) ...
      - ( ( Rei * 1i * k * diag(Wi) + 2.0 * k.^2 * Ii ) * diag(ri.^2) - 3.0 * Ii ) * diag(ri) * D1i(2:end,:) ...
      + ( ( Rei * 1i * k * diag(Wi) * diag(ri.^2) ) * ( k.^2 * diag(ri.^2) + Ii ) + k.^4 * diag(ri.^4) +  2.0 * k.^2 * diag(ri.^2) - 3.0 * Ii ) * D0i;

 Ao = diag(ro.^4) * D4o(1:end-1,:) ...
      + 2.0 * diag(ro.^3) * D3o(1:end-1,:) ...
      - ( ( Reo * 1i * k * diag(Wo) + 2.0 * k.^2 * Io ) * diag(ro.^2) + 3.0 * Io ) * diag(ro.^2) * D2o(1:end-1,:) ...
      - ( ( Reo * 1i * k * diag(Wo) + 2.0 * k.^2 * Io ) * diag(ro.^2) - 3.0 * Io ) * diag(ro) * D1o(1:end-1,:) ...
      + ( ( Reo * 1i * k * diag(Wo) * diag(ro.^2) ) * ( k.^2 * diag(ro.^2) + Io ) + k.^4 * diag(ro.^4) +  2.0 * k.^2 * diag(ro.^2) - 3.0 * Io ) * D0o;

 Bi = -Rei * 1i * k * diag(ri.^2) * ( diag(ri.^2) * D2i(2:end,:) + diag(ri) * D1i(2:end,:) - k.^2 * diag(ri.^2) *  D0i + D0i );
 Bo = -Reo * 1i * k * diag(ro.^2) * ( diag(ro.^2) * D2o(1:end-1,:) + diag(ro) * D1o(1:end-1,:) - k.^2 * diag(ro.^2) *  D0o + D0o );
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% interface BC (5.27b-3): %%%
 Abc3 = zeros(1, Ni+No+1);
 Abc3(1:No+1) -= m * Rei * Wi(1.0) * ( D3o(end,:) + 2.0 * D2o(end,:) - (3.0 * k^2 + 1) * D1o(end,:) - (m - 1.0) * (k^2 - 1.0) * D0o(end,:) / m ) - k * 1i * J *(k^2 - 1.0) * D0o(end,:);
 Abc3(No+1:end) +=   Rei * Wi(1.0) * ( D3i(1,:) + 2.0 * D2i(1,:) - (3.0 * k^2 + 1) * D1i(1,:) );
 Bbc3 = zeros(1, Ni+No+1);
 Bbc3(1:No+1) -= m * Rei * ( D3o(end,:) + 2.0 * D2o(end,:) - (3.0 * k^2 + 1) * D1o(end,:) - (m - 1.0) * (k^2 - 1.0) * D0o(end,:) / m );
 Bbc3(No+1:end) +=   Rei * ( D3i(1,:) + 2.0 * D2i(1,:) - (3.0 * k^2 + 1) * D1i(1,:) );
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 A = [[Ao Oi];Abc3;[Oo Ai]];
 B = [[Bo Oi];Bbc3;[Oo Bi]];
 ii = size(Ao,1) + 1;
 %%% interface BC (5.27b-2): %%% ( using DOF ii+1 )
 D12i = D1i + D2i;
 D12o = D1o + D2o;
 Abc2 = zeros(1, Ni+No+1);
 Abc2(1:No+1) -= m * D12o(end,:);
 Abc2(No+1:end) +=   D12i(1,:);
 Bbc2 = zeros(1, Ni+No+1);
 for p = 1:size(A,2)
  A(ii+1,p) = Abc2(p);
  B(ii+1,p) = Bbc2(p);
 end
 A(ii+1, ii) += (k^2-1.0)*(1.0-m);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% interface BC (5.27b-1): %%% ( using DOF ii-1 )
 Abc1 = zeros(1, Ni+No+1);
 Abc1(1:No+1) -= Wat1 * D1o(end,:);
 Abc1(No+1:end) += Wat1 * D1i(1,:);
 Bbc1 = zeros(1, Ni+No+1);
 Bbc1(1:No+1) -= D1o(end,:);
 Bbc1(No+1:end) += D1i(1,:);
 for p = 1:size(A,2)
  A(ii-1,p) = Abc1(p);
  B(ii-1,p) = Bbc1(p);
 end
 A(ii-1, ii) -= (m - 1.0) * dW2at1;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 [eve,eva] = eig(A,B); % solving A u = B c u
 c = diag(eva);
 omega = c * k;

 %%% filtering spurious eigenvalues: %%%
 tau = imag(omega);
 tau = tau .* isfinite(tau);
 tau = tau .* (tau < MAXeva);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 [taux(ik), ip] = max(tau);
 if (0)%(k<5) % taking second largest
  tau(ip) = -1000000;
  [taux(ik), ip] = max(tau);
 end
 GR = [GR; k taux(ik)];
 %%% plotting eigenvectors: %%%
% for ip = 1:size( c )
% u = eve(:,ip);
% figure(1)
% plot(r, u,'bk'); hold on; plot(r, W, 'r')
% ylabel('w'); xlabel('r');legend('perturbation','base flow')
% pause(1)
% end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%% reproducing Fig. 2 of paper I %%%
figure(2);hold on
plot(kk, taux,'x');xlabel('k');ylabel('Im(\omega)');
title('growth rate');
csvwrite('growRate.csv',GR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (Fall == 1)
 return
end
 [tauMax, ikMax] = max(taux);
 kMax = kk(ikMax);
%%% reproducing Fig. 2 of paper II %%%
 Fig2sig = [Fig2sig; [Re tauMax]];
 Fig2k = [Fig2k; [Re kMax]];
 figure(1);hold on
 plot(Re,tauMax,'dbk')
 xlabel('Re');ylabel('\sigma_{max}');
 figure(2);hold on
 plot(Re,kMax,'^bk');
 xlabel('Re');ylabel('k_{max}');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
csvwrite('sigFig2.csv',Fig2sig)
csvwrite('kFig2.csv',Fig2k)
