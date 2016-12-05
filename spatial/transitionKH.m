% absolute/convective transition of Kelvin-Helmholtz instability
Rt = 1 / 3;
% R < 1 / 3: Im(omega_0) = 0
% R > 1 / 3: None of the omega_0 branches ever crosses the real axis
RR = [-7:0.001:Rt];%RR = [2000:100:10000];
%RR = [-5 -1];
%RR = [0.01:0.01:Rt+0.1];
for l = 1:length(RR)
 R = RR(l);
 k0p = 2 * ( ( 3 + 1/R ) + sqrt( ( 1/R )*( 1/R - 3 ) ) ) / 9;
 k0m = 2 * ( ( 3 + 1/R ) - sqrt( ( 1/R )*( 1/R - 3 ) ) ) / 9;
 figure(1)
 plot(real(k0m),imag(k0m),'.b',real(k0p),imag(k0p),'+bk'); hold on
 xlabel('Re(k_0)');ylabel('Im(k_0)');legend('k_0^-','k_0^+')
 k0 = k0m;
 om0p = k0 + sqrt(R * ( k0^3 - k0^2 ) );
 om0m = k0 - sqrt(R * ( k0^3 - k0^2 ) );
 om0i(l) = abs(imag(om0m)) + abs(imag(om0p));
 figure(2)
 plot(real(om0m),imag(om0m),'.g',real(om0p),imag(om0p),'.r'); hold on
 k0 = k0p;
 om0p = k0 + sqrt(R * ( k0^3 - k0^2 ) );
 om0m = k0 - sqrt(R * ( k0^3 - k0^2 ) );
 figure(2)
 plot(real(om0m),imag(om0m),'+g',real(om0p),imag(om0p),'+r'); hold on
 xlabel('Re(\omega_0)');ylabel('Im(\omega_0)');legend('\omega_0^-','\omega_0^+')
 om0 = om0m;
 om0r(l) = real(om0);
 om0i(l) += abs(imag(om0m)) + abs(imag(om0p));
end
%figure
%plot(RR, om0r, RR, om0i, 'r')
%xlabel('R');legend('Re(\omega_0)','||Im(\omega_0)||')
%max(om0i)
return

kr=[-2:0.1:2];
ki=[-0.9:0.1:0.9];
[kkr,kki] = meshgrid(kr,ki);
k=kkr+i*kki;
omegaP = k + sqrt(R * ( k.^3 - k.^2 ) );
omegaM = k - sqrt(R * ( k.^3 - k.^2 ) );
omeg = omegaM;
%figure;
%surf(kkr,kki,real(omeg))
%xlabel('Re(k)'); ylabel('Im(k)');
%figure;
%surf(kkr,kki,imag(omeg))
%xlabel('Re(k)'); ylabel('Im(k)');
%return

pt = k0m;

k = pt + [-2:0.01:2];
omegaP = k + sqrt(R * ( k.^3 - k.^2 ) );
omegaM = k - sqrt(R * ( k.^3 - k.^2 ) );
omeg = omegaM;
figure% plot how omeg varies with kr
plot(real(k), real(omeg), real(k), imag(omeg), 'r')
hold on;ptr = real(pt) * ones(length(k),1);plot(ptr,omeg,'bk')
xlabel('k');legend('Re(\omega)','Im(\omega)');grid on

k = pt + 1i * [-2:0.01:2];
omegaP = k + sqrt(R * ( k.^3 - k.^2 ) );
omegaM = k - sqrt(R * ( k.^3 - k.^2 ) );
omeg = omegaM;
figure% plot how omeg varies with ki
plot(imag(k), real(omeg), imag(k), imag(omeg), 'r')
hold on;pti = imag(pt) * ones(length(k),1);plot(pti,omeg,'bk')
xlabel('-i k');legend('Re(\omega)','Im(\omega)');grid on
