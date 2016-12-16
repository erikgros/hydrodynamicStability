% comparing results with references
Fall = 1
Ni = 25; % number of points inner region
thres = 0.3; % threshold above which eigenvalues are not considered

if (Fall == 1)
 n = 0;
 a = 1.43;
 m = 0.5;
 zeta = 1;
 J = 0;
 ReiRei = [26.42];
 kk = [0.01:0.5:25];
elseif (Fall == 2)
 eta = 0.8;
 n = 0;
 a = 1.0 / 0.8; % outer radius (inner radius = 1)
 m = 0.1; % viscosity ratio
 zeta = 1;
 J = 0.8*1000.0; % surface tension parameter
 ReiRei = 0.8 * [10:20:500];
 kk = [0.01:0.01:5];
elseif (Fall == 8)
% Fig 8 of paper II:
 eta = 0.95;
 n = 0;
 a = 1.0 / eta; % outer radius (inner radius = 1)
 m = 10.0; % viscosity ratio
 zeta = 1.0;
 J = eta * 100000.0; % surface tension parameter
 ReiRei = 0.95 * [10:500:5000];
 kk = [0.01:0.05:5];
end

sigF = [];
kF = [];
for iRR = 1:length(ReiRei)
 Rei = ReiRei(iRR);
 Re = Rei / eta; % for paper II data

 [taux, ph] = sysI(n, a, m, zeta, J, Rei, Ni, kk);
% taux = sysII(a, m, J, Rei, Ni, kk, thres);

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
 figure(2);hold on
 plot(Re,kMax,'^bk');
 xlabel('Re');ylabel('k_{max}');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
csvwrite(['sigFig' num2str(Fall) '.csv'], sigF)
csvwrite(['kFig' num2str(Fall) '.csv'], kF)
