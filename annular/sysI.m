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
[ri, DMi] = chebdif(Ni+2, 2); % using an additional point at each boundary
ri = 0.5 * (ri + 1.0); % rescaling interval
for n=1:2
 DMi(:,:,n) *= (2.0/1.0)^n; % rescaling derivative
end
D1i=DMi(:,:,1);D2i=DMi(:,:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Matrices (outer region): %%%
[ro, DMo] = chebdif(No+2, 2); % using an additional point at each boundary
DMo = DMo(2:end,2:end,:); % u(a) = v(a) = w(a) = 0
ro = ro(2:end);
ro = Lo * 0.5 * (ro + 1.0) + 1.0; % rescaling interval
for n=1:2
 DMo(:,:,n) *= (2.0/Lo)^n; % rescaling derivative
end
D1o=DMo(:,:,1);D2o=DMo(:,:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = [ro;ri(2:end)];
D0 = eye(Ni+No+3);
Oi = zeros(No+1,Ni+2);
Oo = zeros(Ni+2,No+1);
D1 = [[D1o Oi]; ...
      [Oo D1i]];
D2 = [[D2o Oi]; ...
      [Oo D2i]];

ii = No + 1; %one = r(ii)

%%% Base flow: %%%
[Wi, dWi] = baseFlowI(ri,a,m);
[Wo, dWo] = baseFlowI(ro,a,m);
[W, dW] = baseFlowI(r,a,m);
[Wat1, dW2at1] = baseFlowI( 1.0, a, m);
dWjump = 2.0*(1.0-m)/(a-(1.0-m)/a);
%%%%%%%%%%%%%%%%%%
return
for ik = 1:length(kk)
 k = kk(ik);
end

