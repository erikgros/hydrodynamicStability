% Solving the Eigenvalue problem for the Taylor-Couette instability

N=60; %number interior points

% assembling matrices:
I=eye(N);
Zero=zeros(N);
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
y = 0.5 * y;
halfMinusY = 0.5 * I - diag(y);
for n=1:4 % rescaling:
 DM(:,:,n) *= (2.0)^n;
end
D1=DM(:,:,1);D2=DM(:,:,2);D3=DM(:,:,3);D4=DM(:,:,4);

T = 3000; % initial value Taylor number
for it=1:5
 % secant method
 if (it == 2)
  T = 3800; % second initial value Taylor number
 elseif (it > 2)
  T += -maxk(it-1) * (TT(it-1) - TT(it-2)) /(maxk(it-1) - maxk(it-2));
 end;

 kk=[2:0.1:4];
 for ik=1:length(kk)
  k=kk(ik);
  M = [[ D2-k^2*I Zero ]; [Zero I ] ];
  L = [[ (D4-2*k^2*D2+k^4*I) T*halfMinusY*k^2*I]; [-1*I D2-k^2*I] ];
  eiv=eig(L,M);
  sorted = unique(real(eiv));
  tx(ik)  = sorted( length(sorted) ); % max eigenvalue
  tx2(ik) = sorted( length(sorted)-1 );
  tx3(ik) = sorted( length(sorted)-2 );
  [ev,D]=eig(L,M); % getting the eigenvectors ev and eigenvalued diagonal matrix D

  % finding eigenvector to max (most unstable) eigenvalue:
  imax=-1;
  for iv=1:2*N
   if (real(D(iv,iv)) == tx(ik))
    imax = iv;
   elseif (real(D(iv,iv)) > tx(ik))
    imax = -13
    quit;
   end;
  end;

  % getting normalized eigenvectors for ur and utheta:
  ur=ev(1:N,imax);
  ur = ur / max(abs(ur));
  utheta=ev(N+1:2*N,imax);
  utheta = utheta / max(abs(utheta));

  % finding uz through continuity:
  uz = -(1/k) * D1 * ur;
 end;

 TT(it) = T;
 [maxk(it) ikMax] = max(tx);

 if (abs(maxk(it)) < 0.01)
  display('Found critical Taylor number:');
  T
  kk(ikMax)
  % plotting eigenvectors:
  figure(1)
  plot(y,ur,'b');
  hold all;
  plot(y,utheta,'g');
  plot(y,uz,'r');
  xlabel("y");
  ylabel("Eigenvectors of largest eigenvalue");
  legend("uz","ut","uz");
  break;
 end
%plot(z,sin(z*pi),'ro')
% plot(kk,tx,"b");plot(kk,tx2,"g");plot(kk,tx3,"r");legend("max","2nd","3d");
% plot(kk,tx);
end;
figure(2);
plot(TT, maxk);
hold on
xlabel("Taylor number");
ylabel("max_k(Re(eigenvalues))");
