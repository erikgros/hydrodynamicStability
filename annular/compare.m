% comparing results with references
Fall = 2
Ni = 10; % number of points inner region
No = Ni;

n = 0;
zeta = 1.0;
E = 0; % !computing energy
if (Fall == 1)
 a = 1.43;
 m = 0.5;
 J = 0;
 ReiRei = [26.42];
 kk = [0.01:0.5:25];
elseif (Fall == 2)
 eta = 0.8;
 E = 1;
 a = 1.0 / eta; % outer radius (inner radius = 1)
 m = 0.1; % viscosity ratio
 J = eta*1000.0; % surface tension parameter
 ReiRei = eta * [10:20:500];
 kk = [0.01:0.01:3];
 if (E == 1 && eta == 0.8)
  kk = [7.00000000000000e-01 6.80000000000000e-01 6.40000000000000e-01 6.00000000000000e-01 5.40000000000000e-01 4.50000000000000e-01 1.00000000000000e-02 1.00000000000000e-02 1.00000000000000e-02 2.63000000000000e+00 2.51000000000000e+00 2.40000000000000e+00 2.31000000000000e+00 2.25000000000000e+00 2.19000000000000e+00 2.14000000000000e+00 2.10000000000000e+00 2.06000000000000e+00 2.02000000000000e+00 1.99000000000000e+00 1.95000000000000e+00 1.92000000000000e+00 1.89000000000000e+00 1.87000000000000e+00 1.84000000000000e+00];
 elseif (eta == 0.7)
  kk=[0.3300000000000000 0.6900000000000001 0.6800000000000000 0.6600000000000000 0.6400000000000000 0.6200000000000000 0.6100000000000000 0.6100000000000000 0.6400000000000000 0.7100000000000001 0.7700000000000000 0.8000000000000000 0.8100000000000001 0.8100000000000001 0.8000000000000000 0.7800000000000000 0.7700000000000000 0.7500000000000000 0.7300000000000000 0.0400000000000000 0.0400000000000000 0.0400000000000000 0.0400000000000000 0.0400000000000000 0.0400000000000000];
 end
 Ni = 20;
 No = 30;
elseif (Fall == 8)
% Fig 8 of paper II:
 eta = 0.95;
 a = 1.0 / eta; % outer radius (inner radius = 1)
 m = 10.0; % viscosity ratio
 J = eta * 100000.0; % surface tension parameter
 ReiRei = eta * [10 80 200 800 2000 3500 5000];
 kk = [0.01:0.05:5];
elseif (Fall == 11)
%Fig 11 of paper II
 eta = 0.2;
 a = 1.0 / eta; % outer radius (inner radius = 1)
 m = 0.1; % viscosity ratio
 J = eta*1000.0; % surface tension parameter
% m = 10.0; % viscosity ratio
% J = eta*100000.0; % surface tension parameter
 ReiRei = eta * [5 50 100 200 500];
 kk = [0.01:0.01:1];
end

sigF = [];
kF = [];
for iRR = 1:length(ReiRei)
 Rei = ReiRei(iRR);
 Re = Rei * a; % for paper II data

 if ( E == 1 )
  k = kk(iRR);
  [taux, ph] = sysI(n, a, m, zeta, J, Rei, Ni, No, k);
 else
  [taux, ph] = sysI(n, a, m, zeta, J, Rei, Ni, No, kk);
% [taux, ph] = sysII(a, m, J, Rei, Ni, No, kk);
 end

 if (Fall == 1)
  GR = [kk' taux'];
 %%% reproducing Fig. 2 of paper I %%%
% figure(3);hold on
% plot(kk, taux,'x');xlabel('k');ylabel('Im(\omega)');
% title('growth rate');
  csvwrite(['growRate' num2str(n) '.csv'],GR)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  return
 end
 [tauMax, ikMax] = max(taux);
 kMax = kk(ikMax);
%%% reproducing Fig. 2 of paper II %%%
 sigF = [sigF; [Re tauMax]];
 kF = [kF; [Re kMax]];
 figure(1);hold on
 plot(Re,tauMax,'dbk')
 xlabel('Re');ylabel('\sigma_{max}');
set(gca,'xscale','log');
set(gca,'yscale','log');
 figure(2);hold on
 plot(Re,kMax,'^bk');
 xlabel('Re');ylabel('k_{max}');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
csvwrite(['sigFig' num2str(Fall) '.csv'], sigF)
csvwrite(['kFig' num2str(Fall) '.csv'], kF)
