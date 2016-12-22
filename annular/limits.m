% comparing limiting cases
Ni = 10; % number of points inner region
n = 0;
J = 1e10;
k = [0.001:0.05:1.001];
%% Rayleigh-Plateau: %%
No = 6;
eta = 0.4
m = 1e-3;
zeta = 1e-4;
Rei = 1e-3;
factor = sqrt(J) / Rei;
omega = sqrt( k .* (k.^2 - 1) .* besseli(1,k) ./ besseli(0,k) );
name = "Rayleigh-Plateau";
%% Chandrasekhar: %%
No = 30;
eta = 0.1
m = 1000;
zeta = 10000;
Rei = 1e-2;
factor = sqrt(J/zeta) / Rei;
omega = sqrt( k .* (k.^2 - 1) .* besselk(1,k) ./ besselk(0,k) );
name = "Chandrasekhar";
%%%%%%%%%%%%%%%%%%%%

a = 1.0 / eta; % outer radius (inner radius = 1)
[taux, ph] = sysI(n, a, m, zeta, J, Rei, Ni, No, k);
 %%% plotting growth rate: %%%
 figure(1);hold on
 plot(k, taux,'+g');xlabel('k');ylabel('Im(\omega)');
 title('growth rate');
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tauA = factor * imag(omega);
plot(k, tauA,'bk')
legend("",name)
xlim([0 1.001])
