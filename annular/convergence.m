% benchmarking results with Table 1 of paper II
Case = 1

% Case 2:
Jstar = 0.0;
eta = 0.7;
m = 0.5;
zeta = 1;
n = 0;
beta = 10.0;
Re = 37.78;

if (Case == 1)
 Jstar = 1000.0;
 eta = 0.9;
 m = 0.05;
 beta = 5.0;
 Re = 500.0;
elseif (Case == 3)
 n = 5;
end

J = eta * Jstar;
a = 1.0 / eta;
k = eta * beta;
Rei = eta * Re;

NiNi = [10 20 40 80 160 320]; % # points inner region
for l = 1:length(NiNi)
 Ni = NiNi(l);
 [tau, c(l) ] = sysI(n, a, m, zeta, J, Rei, Ni, k);
 c(l)
end
err = abs(c .- c(end));
semilogy(NiNi, err, 'g')
hold on
xlabel('N_i');ylabel('error');
title('convergence')
