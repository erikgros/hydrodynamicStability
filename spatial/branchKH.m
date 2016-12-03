% finding spatial branches of the Kelvin-Helmholtz dispersion relation
%RR = [0.1:0.1:1.1];%
RR=[0.08]; % slip ratio
for l=1:length(RR)
 R = RR(l);
 om = [0.001:0.05:1.051];
 %%% Gaster transform: %%%
 for j = 1:length(om)
  k = om(j);
  pp1 = [1 -2*k k*k*(R+1)-R*k^3];
  r1 = roots( pp1 ); % roots of quadratic polynomial
  om1(j) = max( imag(r1) ); % = imag(k + sqrt(R * ( k.^3 - k.^2 ) ))
 end
 figure(2);
 plot(om, om1,'r'); hold all;
 %%%%%%%%%%%%%%%%%%%%%%%%%
 for j = 1:length(om)
  ppp = [R -(R+1) 2*om(j) -(om(j))^2];
  r(j,:) = roots(ppp); % roots of cubic polynomial
 end

 for k=1:3 % plotting the 3 roots/branches:
  ki = imag(r(:,k));
  kr = real(r(:,k));
  figure(1);
  plot(kr, ki,'x'); hold all;
  xlabel('k_r'); ylabel('k_i');
  figure(2);
  plot(om, -ki,'x'); hold all;
  xlabel('\omega'); ylabel('-k_i');
  figure(3);
  plot(om, kr,'x'); hold all;
  xlabel('\omega'); ylabel('k_r');
 end
end;
