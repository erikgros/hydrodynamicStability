% OPTIMAL EXCITATION: initial condition maximizing the gain (G) at a given time
Tmax = 5;

% Charru sect.1.5.1:
%epsi = 0.1;
%L = [ -epsi, 1; 0 -epsi*2];
%a = 3 + sqrt( 1-8*epsi^2 ); a = a*0.25 / ( 1 + epsi^2 )
% Charru sect.1.5.2:
%epsi1 = 0.1;epsi2 = epsi1;
%L = [ epsi1, 1; 0 -epsi2];
% exercice on transient growth:
for Re = 1:10:101%Re = 1
LS = (1/Re)*[-3, 1; 1,-5];
%L = LS; % Question 1
LOS = (1/Re)*[0, -1; 4,-5];
%L = LOS; % Question 2
L = [[LOS zeros(2,2)];[ones(2) LS]]; % Question 3

[eve,eva] = eig(L);
[evea,evaa] = eig(L'); % adjoint matrix
[evamax, imax] = max(diag(evaa))
[evamin, imin] = min(diag(eva))
%prodEVE = eve(:,1)' * eve(:,2)
M = (eve^-1)' * eve^-1; % Question 7

% 2 initial conditions:
x0maxLT = evea(:,imax); %  maximum growth at long times
%orthogonal = x0maxLT' * eve(:,imin)
%y0 = [1; a];
%x0maxLT = y0(1) * eve(:,1) + y0(2) * eve(:,2);
%x0maxLT = x0maxLT / sqrt(x0maxLT(1)^2+x0maxLT(2)^2);

T = Tmax * Re;
tt = [0:0.01*T:T];
for it = 1:length(tt)
 t = tt(it);
 expLt = expm( L * t );
 x = expLt * x0maxLT;
 G(it) = sqrt(x' * x);
 % finding optimal perturbation by computing svd of expLt
 %% long story:
 SA = expLt' * expLt;
 [sve,ssva] = eig(SA); % ssva are squared singular values of expLt
 [svamax, simax] = max(diag(sqrt(ssva)));
 x0max = sve(:,simax);
 [U, S, V] = svd(expLt); svamax = max(diag(S));
 x0max = V(:,1); % sorted in decreasing order
 x = expLt * x0max;
 Gopt(it) = sqrt(x' * x);
 % Question 7:
 GMopt(it) = sqrt(x' * M * x);
% [U, S, V] = svd((eve^-1) * expLt); svamax = max(diag(S));
% GMopt(it) = svamax;
 %% short story:
 Gopt(it) = norm(expLt);
 GMopt(it) = norm((eve^-1) * expLt); % Question 7
end
maxGain = max( G );
% Question 6:
%figure(1)
%plot( tt/Re, log(Gopt/Re)),shg
%hold on;
%ylabel('log(Gopt/Re)');xlabel('t/Re');
%figure(2)
%plot(Re, maxGain/Re)
%hold on;
%ylabel('maxG/Re');xlabel('Re');
figure
plot( tt, log(Gopt));
hold on;
ylabel('log(G)');xlabel('t');
plot( tt, log(G), 'r');
hold on;
plot(tt, evamax*tt,'g'); legend('G_{opt}(t)','G_{opt}(\infty)','maxeva');
end % end Re loop
