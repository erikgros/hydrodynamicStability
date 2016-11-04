% Solving the Eigenvalue problem for the Taylor-Couette instability

N=60; %number interior points
dz=1/(N+1);
y=[-0.5+dz:dz:0.5-dz];

% assembling matrices:
I=eye(N);
halfMinusY = 0.5 * I - diag(y);
Zero=zeros(N);
D1=Zero;
for ii=1:N-1
 D1(ii,ii+1) =  1;
 D1(ii+1,ii) = -1;
end
D1=D1/dz;
D2=-2*eye(N);
for ii=1:N-1
 D2(ii,ii+1) = 1;
 D2(ii+1,ii) = 1;
end
D2=D2/dz^2;
D4=6*eye(N);
for ii=1:N-1
 D4(ii,ii+1) = -4;
 D4(ii+1,ii) = -4;
 if ( ii <= N-2)
  D4(ii,ii+2) = 1;
  D4(ii+2,ii) = 1;
 end
end

% no-slip BC + continuity:
D4(1,1) = 7;
D4(N,N) = 7;

D4=D4/dz^4;

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
 maxk(it) = max(tx);

 if (abs(maxk(it)) < 0.01)
  display('Found critical Taylor number:');
  T
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
xlabel("Taylor number");
ylabel("max_k(Re(eigenvalues))");
