% finding spatial branches of the Kelvin-Helmholtz dispersion relation
omr = [0.001:0.05:1.051];
R = 0.1; % slip ratio
tauMax = 2*sqrt(R/3) / 3; % max growth rate for any real k
gamma = [0.0:0.2*tauMax:1.4*tauMax];
% using fact that there cannot be a real k with Im(omega) > tauMax to distinguish branches:
% the branch that cuts the real axis at gamma = tauMax is the one corresponding to Im(omega) > 0
for l=1:length(gamma)
 om = omr + i*gamma(length(gamma)-l+1);

 %%% Gaster transform: %%%
 for j = 1:length(om)
  k = omr(j);
  pp1 = [1 -2*k k*k*(R+1)-R*k^3];
  r1 = roots( pp1 ); % roots of quadratic polynomial
  om1(j) = max( imag(r1) ); % = imag(k + sqrt(R * ( k.^3 - k.^2 ) ))
 end
% figure(2); plot(omr, om1,'r'); hold all;
 %%%%%%%%%%%%%%%%%%%%%%%%%
 for j = 1:length(om)
  ppp = [R -(R+1) 2*om(j) -(om(j))^2];
  r(j,:) = roots(ppp); % roots of cubic polynomial
 end
 figure
 for k=1:3 % plotting the 3 roots/branches:
  ki = imag(r(:,k));
  kr = real(r(:,k));
%  figure(1);
  plot(kr, ki,'x'); hold all;
  xlabel('k_r'); ylabel('k_i');
  grid on;
%  figure(2);
%  plot(omr, -ki,'x'); hold all;
%  xlabel('\omega'); ylabel('-k_i');
%  figure(3);
%  plot(omr, kr,'x'); hold all;
%  xlabel('\omega'); ylabel('k_r');
 end
end;
