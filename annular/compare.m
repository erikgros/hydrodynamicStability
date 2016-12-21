% comparing results with references
Fall = 11
Ni = 10; % number of points inner region

n = 0;
zeta = 1.0;
if (Fall == 1)
 eta = 100000000
 a = 1.43;
 m = 0.5;
 J = 0;
 ReiRei = [26.42];
 kk = [0.01:0.5:25];
elseif (Fall == 2)
 eta = 0.8;
 a = 1.0 / eta; % outer radius (inner radius = 1)
 m = 0.1; % viscosity ratio
 J = eta*1000.0; % surface tension parameter
 ReiRei = eta * [10:20:500];
 kk = [0.01:0.01:5];
elseif (Fall == 8)
% Fig 8 of paper II:
 eta = 0.95;
 a = 1.0 / eta; % outer radius (inner radius = 1)
 m = 10.0; % viscosity ratio
 J = eta * 100000.0; % surface tension parameter
 ReiRei = eta * [10 20 40 90 180 360 720 1440];
%               [10 80 200 800 2000 3500 5000];
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
 Re = Rei / eta; % for paper II data

 [taux, ph] = sysI(n, a, m, zeta, J, Rei, Ni, kk);
% [taux, ph] = sysII(a, m, J, Rei, Ni, kk);

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
