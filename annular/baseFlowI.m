function [W, dWdr] = baseFlowI( r, a, m )
% supply r, a and density ratio m
eta=1/a;
 A = 1.0-(1.0-m)/a^2;
 s = (r < 1.0);
 W = s.*(1.0-eta^2+m*(1.0-r.^2)/a^2)/A + (!s).*(1.0-(r/a).^2)/A;
 dWdr = -2.0*(s*m+!s).*r/(a*A);
end
