% solving Rayleigh's equation to study the stability of a tanh velocity profile
% in inviscid flow between two parallel plates
N = 80;
L = 5; % half of distanve between plates
dy = 2*L/N;
y = [-L:dy:L]';
I = eye(N-1);
dM = -2*ones(N-1,1);
D2v = diag(dM);
for ip=2:N-2;
  D2v(ip,ip-1)=1;
  D2v(ip,ip+1)=1;
end;
D2v(1,2)=1;
D2v(N-1,N-2)=1;
D2v=D2v/dy^2;

Up = 1 + tanh(y(2:N)); % tanh profile
U2p = -2*sinh(y(2:N))./cosh(y(2:N)).^3; % second derivative of tanh

kk=[0.02:0.02:1.02];%kk=0.3;
for p=1:length(kk)
  k = kk(p);
  A = ( diag(Up) * (D2v - k^2 * I) - diag(U2p));
  B = D2v - k^2 * I;
  [evec,eva] = eig(A,B);
  l = diag(eva) * k; % eva is c (= omega / k)
  [taux(p), ip] = max(imag(l));
  freq(p) = real( l(ip) );
  figure(1);
  %plot(real(l),imag(l),'x');
  %xlabel('real(\omega)');
  %ylabel('imag(\omega)');

  % eigenvectors:
  phi = evec(:,ip);
  phi = phi / max(abs(phi));
  plot(y(2:N),phi,'r');
  xlabel("y");
  ylabel("Most unstable perturbation");
  legend('\psi(y)');
end;
figure(2);hold all;
plot(kk, taux(1:length(kk)),'o');
xlabel('k');
ylabel('imag(\omega)');
figure(3);hold all;
plot(kk,freq(1:length(kk)),'o');
xlabel('k');
ylabel('re(\omega)');
