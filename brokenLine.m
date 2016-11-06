% last question of the exercice about broken line profile
U1 = 1;
U2 = 0;
dU = U1 - U2;
Ubar = 0.5 * ( U1 + U2 );

% finding max kDelta for instability:
x = 1;
res = x - 1 - exp(-x);
while( abs(res) > 0.0000001 )
 dres = 1 + exp(-x);
 x = x - res / dres;
 res = x - 1 - exp(-x);
end
x

% computing kDelta for max growth rate:
kDeltaMax = 1;
res = kDeltaMax - 1 + exp(-2*kDeltaMax);
while( abs(res) > 0.0000001 )
 dres = 1 - exp(-2*kDeltaMax);
 kDeltaMax = kDeltaMax - res / dres;
 res = kDeltaMax - 1 + exp(-2*kDeltaMax);
end
kDeltaMax

kDelta = [x-1.0 : 0.1 : x+10];
discrim = ((kDelta .- 1).^2 .- exp(-2 .* kDelta)) .* (dU ./ kDelta).^2;
plus = kDelta .* Ubar .+ 0.5 .* kDelta .* sqrt(discrim);
figure(1)
plot(kDelta, real(plus))
hold on;
minus = kDelta .* Ubar .- 0.5 .* kDelta .* sqrt(discrim);
plot(kDelta, real(minus), 'r')
xlabel('k \delta');
ylabel('Real part');
legend('(\omega \delta)_+','(\omega \delta)_-')
figure(2)
plot(kDelta, imag(plus))
xlabel('k \delta');
ylabel('growth rate: Im((\omega \delta)_+)');
