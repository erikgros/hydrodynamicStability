% comparing limiting cases
Ni = 10; % number of points inner region
% eta = 0.1;
eta = 0.8
n = 0;
a = 1.0 / eta; % outer radius (inner radius = 1)
large = 1e10
% m = mid; % viscosity ratio
% zeta = mid;
m=1e-3;
zeta=1e-4;
J=large;
ReiRei = [1e-3];
k = [0.01:0.01:1];

for iRR = 1:length(ReiRei)
 Rei = ReiRei(iRR);
% Re = Rei / eta; % for paper II data

 [taux, ph] = sysI(n, a, m, zeta, J, Rei, Ni, k);

 %%% plotting growth rate: %%%
 figure(1);hold on
 plot(k, taux,'r');xlabel('k');ylabel('Im(\omega)');
 title('growth rate');
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 [tauMax, ikMax] = max(taux);
 kMax = k(ikMax);
end
omega = sqrt( k .* (k.^2 - 1) .* besseli(1,k) ./ besseli(0,k) );
tauRP = imag(omega);
factor = tauMax / max( tauRP);
plot(k,factor*tauRP,'bk')
