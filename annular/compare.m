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
 n = 0;
 a = 1.0 / 0.8; % outer radius (inner radius = 1)
 m = 0.1; % viscosity ratio
 zeta = 1;
 J = 0.8*1000.0; % surface tension parameter
 ReiRei = 0.8 * [10:20:500];
 kk = [0.01:0.01:5];
end

Fig2sig = [];
Fig2k = [];
for iRR = 1:length(ReiRei)
 Rei = ReiRei(iRR);
 Re = Rei / 0.8; % for paper II data

 taux = sysI(n, a, m, zeta, J, Rei, Ni, kk);
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
 Fig2sig = [Fig2sig; [Re tauMax]];
 Fig2k = [Fig2k; [Re kMax]];
 figure(1);hold on
 plot(Re,tauMax,'dbk')
 xlabel('Re');ylabel('\sigma_{max}');
 figure(2);hold on
 plot(Re,kMax,'^bk');
 xlabel('Re');ylabel('k_{max}');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
csvwrite('sigFig2.csv',Fig2sig)
csvwrite('kFig2.csv',Fig2k)
