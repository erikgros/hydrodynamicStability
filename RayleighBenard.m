% Solving the Eigenvalue problem for the Rayleigh-Benard instability
Pr=1;

N=90; %number interior points
dz=1/(N+1);
z=[dz:dz:N*dz]; % z-coordinate

% assembling matrices:
I=eye(N);
Z=zeros(N);
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
%D4(1,1)=7;
%D4(N,N)=7;

% stress free BC + continuity:
D4(1,1)=5;
D4(N,N)=5;

D4=D4/dz^4;

%% Eigenvalue problem:  s M phi = L phi

%    [ D^2-k^2 I ; Z ] [uz']   [Pr(D^2-k^2 I)^2 ; -Pr k^2 I ] [uz']
%  s |               | |   | = |                            | |   |
%    [    Z      ; I ] [th']   [      Ra I      ; D^2-k^2 I ] [th']

Ra = 1500; % initial value Rayleigh number
for it=1:6
 % secant method
 if (it == 2)
  Ra = 1800; % second initial value Taylor number
 elseif (it > 2)
  Ra += -maxk(it-1) * (TT(it-1) - TT(it-2)) /(maxk(it-1) - maxk(it-2));
 end;

  kk=[1.5:0.1:4]; % all k-values
  for ik=1:length(kk)
    k=kk(ik);
    M=[[D2-k^2*I Z];[Z I ]];
    L=[[Pr*(D4-2*k^2*D2+k^4*I) -Pr*k^2*I];[Ra*I D2-k^2*I]];
    eiv=eig(L,M);
    sorted = unique(real(eiv));
    tx(ik)  = sorted( length(sorted) ); % max eigenvalue
    tx2(ik) = sorted( length(sorted)-1 );
    tx3(ik) = sorted( length(sorted)-2 );

%    hold all;
%    plot(real(eiv),imag(eiv),'ro');

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

    % getting normalized eigenvectors for uz and theta:
    uz=ev(1:N,imax);
    uz = uz / max(abs(uz));
    theta=ev(N+1:2*N,imax);
    theta = theta / max(abs(theta));
 end;
 TT(it) = Ra;
 maxk(it) = max(tx);

 if (abs(maxk(it)) < 0.01)
  display('Found critical Rayleigh number:');
  Ra
  hold all;
  plot(kk,tx,"b");plot(kk,tx2,"g");plot(kk,tx3,"r");legend("max","2nd","3d");
  xlabel('k');ylabel('Re(s)');
  % plotting eigenvectors:
%  figure(1)
%  hold all;
%  plot(z,theta,'g');
%  plot(z,uz,'r');
%  xlabel("z");
%  ylabel("Eigenvectors of largest eigenvalue");
%  legend("theta","uz");
%plot(z,sin(z*pi),'ro')
  break;
 end
end;
%plot(TT, maxk);
%xlabel("Rayleigh number");
%ylabel("max_k(Re(eigenvalues))");
