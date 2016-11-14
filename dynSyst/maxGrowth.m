% OPTIMAL EXCITATION: initial condition maximizing the "energy" at a given time
epsi = 0.1;
L = [ -epsi, 1; 0 -epsi*2]; % Charru sect.1.5.1
%epsi1 = 0.1;epsi2 = epsi1;
%L = [ epsi1, 1; 0 -epsi2]; % Charru sect.1.5.2
%x0 = [1, 0]'; % direction of first eigenvector

[eve,eva] = eig(L);
[evea,evaa] = eig(L'); % adjoint matrix
[evamax, imax] = max(diag(evaa))
[evamin, imin] = min(diag(eva))

% 2 initial conditions:
x0maxLT = evea(:,imax) %  maximum growth at long times
orthogonal = x0maxLT' * eve(:,imin)

tt = [0:0.2:60];
for it = 1:length(tt)
 t = tt(it);
 expLt = expm( L * t );
 SA = expLt' * expLt; % eigenvalues of SA are square of singular values of expLt
 [seve,seva] = eig(SA);
 [sevamax, simax] = max(diag(seva));
 x0max = seve(:,simax);
 x = expLt * x0maxLT;
 E(it) = x' * x;
 x = expLt * x0max;
 Emax(it) = x' * x;
end
plot( tt, log(E), 'r')
hold on;
plot( tt, log(Emax))
ylabel('E');xlabel('t');
legend('maxLT','max');
