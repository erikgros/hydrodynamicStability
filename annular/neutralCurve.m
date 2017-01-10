% Fig 7 paper II

Ni = 15; % number of points inner region
No = 8;

eta = 0.99;
n = 0;
a = 1.0 / eta; % outer radius (inner radius = 1)
J = eta * 100000.0; % surface tension parameter
m = 10.0; % viscosity ratio
zeta = 1.0;
ReiRei = eta * [10 80 200 800 2000 3500 5000];

for iRR = 1:length(ReiRei)
 Rei = ReiRei(iRR);

 k = 0.9;
% it = 0;
 sigma = 1;
 while ( abs( sigma ) > 0.00001 )
  % secant method
%  it += 1;
%  if (it == 2)
%   k = 8; % second initial value
%  elseif (it > 2)
%   k += -sig(it-1) * (kk(it-1) - kk(it-2)) / (sig(it-1) - sig(it-2));
%  end;
  k += 0.001;

  [sigma, dummy] = sysI(n, a, m, zeta, J, Rei, Ni, No, k);
%  kk( it ) = k;
%  sig( it ) = sigma;
 end
 kkk( iRR ) = k;
end

Re = ReiRei / eta; % for paper II data
figure(1);hold on
semilogy(kkk, Re)
xlabel('k');ylabel('Re');
