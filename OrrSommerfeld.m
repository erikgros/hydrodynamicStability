% solving the Orr Sommerfeld equation to study the stability of 2d Poiseuille flow
Re = 5772.22;
N = 60;
h = 1.0; % half of distance between plates

I = eye(N);
Z = zeros(N);
addpath('./Chebyshev')
[y, DM] = chebdif(N+4, 4);
% 0th and first derivative vanish at boundary points:
for n=1:4
 for o=3:N+2
  for p=3:N+2
   DM(p,o,n) -= DM(2,o,1) .* DM(p,2,n) / DM(2,2,1);
   DM(p,o,n) -= DM(end-1,o,1) .* DM(p,end-1,n) / DM(end-1,end-1,1);
  end
 end
end
DM = DM(3:end-2,3:end-2,:); % removing BC points
y = y(3:end-2);
% rescaling interval
y = h * y;
halfMinusY = 0.5 * I - diag(y);
for n=1:4 % rescaling:
 DM(:,:,n) *= (1/h)^n;
end
D1=DM(:,:,1);D2=DM(:,:,2);D3=DM(:,:,3);D4=DM(:,:,4);

Up = 1.0 *( h^2 .- y.^2) / (h^2); % 2d Poiseuille profile
U2p = -2.0 * ones(N,1) / (h^2); % second derivative

kk=[0.01:0.01:1.8];%kk=[1.02056];
tauMax = -1000000000;
for Re = [1772, 5772, 9772]
Re
for p=1:length(kk)
  k = kk(p);
  A = ( diag(Up) * (D2 - k^2 * I) - diag(U2p) - (D4 - 2*k^2*D2 + k^4*I) / (1i*k*Re) );
  B = D2 - k^2 * I;
  [eve,eva] = eig(A,B);
  c = diag(eva);
  omega = c * k;
  [taux(p), ip] = max(imag(omega));
  freq(p) = real( omega(ip) );
  if (taux(p) > tauMax)
   tauMax = taux(p);
   kMax = k;
   omegaMax = omega;
   eveMax = eve(:,ip);
  end
end;
kMax

% plotting:
figure(1)
hold on;
plot(imag(omegaMax),-real(omegaMax),'+');%axis([-1 0 -1 0]); % like Wikipedia
ylabel('-Re(\omega)');xlabel('Im(\omega)');
%plot(real(omegaMax),imag(omegaMax),'o');xlabel('Re(\omega)');ylabel('Im(\omega)');
title("critical conditions")

figure(2);hold all;
plot(kk, taux,'+');xlabel('k');ylabel('Im(\omega)');
end
legend('Re = 1772','     5772','     9772');

%figure;hold all;
%plot(kk,freq);xlabel('k');ylabel('Re(\omega)');

%figure(3)
%phi = eveMax;
%phi = phi / max(abs(phi));
%plot(y,phi,'r');
%xlabel("y");
%ylabel("Streamfunction");
%legend("Most unstable perturbation");
