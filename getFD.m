function [z, DM] = getFD(N)
% finite difference matrices:
dz = 2.0 / (N-1);
z=[1.0:-dz:-1.0]';

D1=zeros(N);
for ii=1:N-1
 D1(ii,ii+1) =  1;
 D1(ii+1,ii) = -1;
end
% using back/forward differences at boundaries:
D1(1,1)=-2;D1(1,2)=2;D1(N,N)=2;D1(N,N-1)=-2;
D1=D1/(2*dz);
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
DM(:,:,1) = -D1;
DM(:,:,2) = D2;
DM(:,:,3) = -D2 * D1;
DM(:,:,4) = D4;

