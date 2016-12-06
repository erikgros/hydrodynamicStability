% solving the Orr Sommerfeld equation to study the stability of 2d Poiseuille flow
Re = 5772.22;
N = 60;
h = 1.0; % half of distance between plates

I = eye(N);
Z = zeros(N);
addpath('./Chebyshev')
[y,DM] = chebdif(N+2,2);
D2 = DM(2:N+1,2:N+1,2); % Dirichlet BCs
%y = x(2:end-1);
[y,D4] = cheb4c(N+2); % using cheb4c to enforce clamped BC
% rescaling interval:
D2(:,:) *= (1/h)^2; % rescaling
D4(:,:) *= (1/h)^4; % rescaling
y = h * y;

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
