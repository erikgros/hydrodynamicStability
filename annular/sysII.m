% axisymetric system II from "Lubricated pipelining: stability of core annular flow"
% only 1 interface BC
a = 1.43; % outer radius (inner radius = 1)
m = 0.5; % viscosity ratio
J = 0.0; % surface tension parameter

Ni = 40; % number interior points
No = 20; % number exterior points

Rei = 26.42;
Reo = m * Rei;
Lo = (a - 1.0);

% assembling matrices:
Ii = eye(Ni);
D0i = [zeros(Ni,1) Ii];
Io = eye(No);
D0o = [Io zeros(No,1)];
Oo = zeros(Ni,No);
Oi = zeros(No,Ni);
addpath('../Chebyshev')
[ri, DMi] = chebdif(Ni+2, 4); % using 2 ghost points (one at each boundary)
DMi = DMi(1:end-1,1:end-1,:); % u1(0) = 0
ri = ri(2:end-1);
ri = 0.5 * (ri + 1.0); % rescaling interval
for n=1:4
 DMi(:,:,n) *= (2.0/1.0)^n; % rescaling derivative
end
D1i=DMi(:,:,1);D2i=DMi(:,:,2);D3i=DMi(:,:,3);D4i=DMi(:,:,4);
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

% base flow:
[Wi, dWi] = baseFlowI(ri,a,m);
[Wo, dWo] = baseFlowI(ro,a,m);
r = [ro;1.0;ri];
kk = [0:0.5:20];
for ik = 1:length(kk)
 k = kk(ik);
 % A, B given by (5.25):
 Ai = diag(ri.^4) * D4i(2:end,:) + 2.0 * diag(ri.^3) * D3i(2:end,:) - ( ( Rei * 1i * k * diag(Wi) + 2.0 * k.^2 * Ii ) * diag(ri.^2) + 3.0 * Ii ) * diag(ri.^2) * D2i(2:end,:) - ( ( Rei * 1i * k * diag(Wi) + 2.0 * k.^2 * Ii ) * diag(ri.^2) - 3.0 * Ii ) * diag(ri) * D1i(2:end,:) + ( ( Rei * 1i * k * diag(Wi) * diag(ri.^2) ) * ( k.^2 * diag(ri.^2) + Ii ) * D0i + k.^4 * diag(ri.^4) * D0i +  2.0 * k.^2 * diag(ri.^2) * D0i- 3.0 * D0i );
 Ao = diag(ro.^4) * D4o(1:end-1,:) + 2.0 * diag(ro.^3) * D3o(1:end-1,:) - ( ( Reo * 1i * k * diag(Wo) + 2.0 * k.^2 * Io ) * diag(ro.^2) + 3.0 * Io ) * diag(ro.^2) * D2o(1:end-1,:) - ( ( Reo * 1i * k * diag(Wo) + 2.0 * k.^2 * Io ) * diag(ro.^2) - 3.0 * Io ) * diag(ro) * D1o(1:end-1,:) + ( ( Reo * 1i * k * diag(Wo) * diag(ro.^2) ) * ( k.^2 * diag(ro.^2) + Io ) * D0o + k.^4 * diag(ro.^4) * D0o +  2.0 * k.^2 * diag(ro.^2) * D0o- 3.0 * D0o );
 Bi = Rei * 1i * k * diag(ri.^2) * ( diag(ri.^2) * D2i(2:end,:) + diag(ri) * D1i(2:end,:) + k.^2 * diag(ri.^2) *  D0i + D0i );
 Bo = Reo * 1i * k * diag(ro.^2) * ( diag(ro.^2) * D2o(1:end-1,:) + diag(ro) * D1o(1:end-1,:) + k.^2 * diag(ro.^2) *  D0o + D0o );
 %%%% interface BC only (5.27b): %%%
 Abc = zeros(1, Ni+No+1);
 Abc(1:No+1) -= m * Rei * Wi(1.0) * ( D3o(end,:) + 2.0 * D2o(end,:) - (3.0 * k^2 + 1) * D1o(end,:) - (m - 1.0) * (k^2 - 1.0) * D0o(end,:) / m ) - k * 1i * J *(k^2 - 1.0) * D0o(end,:);
 Abc(No+1:end) +=   Rei * Wi(1.0) * ( D3i(1,:) + 2.0 * D2i(1,:) - (3.0 * k^2 + 1) * D1i(1,:) );
 Bbc = zeros(1, Ni+No+1);
 Bbc(1:No+1) -= m * Rei * ( D3o(end,:) + 2.0 * D2o(end,:) - (3.0 * k^2 + 1) * D1o(end,:) - (m - 1.0) * (k^2 - 1.0) * D0o(end,:) / m );
 Bbc(No+1:end) +=   Rei * ( D3i(1,:) + 2.0 * D2i(1,:) - (3.0 * k^2 + 1) * D1i(1,:) );
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 A = [[Ao Oi];Abc;[Oo Ai]];
 B = [[Bo Oi];Bbc;[Oo Bi]];
 %A u = B c u
 [eve,eva] = eig(A,B);
 c = diag(eva);
 omega = c * k;
 [taux(ik), ip] = max(imag(omega));
 u = eve(:,ip);
 plot(r, u)
end
figure
plot(kk, taux);xlabel('k');ylabel('Im(\omega)');
title('growth rate');
