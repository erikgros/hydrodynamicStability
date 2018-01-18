% growth rate Vs. wave number
% experiment 3 (Fig 14) paper I:
n = 0;
zeta = 1.0;
a = 1.42;
J = 2102 / a;
m = 0.0532;
Re1 = 69.8;
k = 0.6;

for f=10:10:30
 Ni = 2*f;
 No = f;
 [taux, ph] = sysI(n, a, m, zeta, J, Re1, Ni, No, k);
 taux % should converge to 0.0783518 (Tab. 1 paper I)
end

%%% plotting growth rate: %%%
k = [0.1:0.1:1];
[taux, ph] = sysI(n, a, m, zeta, J, Re1, Ni, No, k);
plot(k, taux,'*-');xlabel('k');ylabel('Im(\omega)');
title('growth rate')
