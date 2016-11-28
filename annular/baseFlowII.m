function [W, dWdr] = baseFlowII( r, eta, m )
% supply r,eta in [0, 1] and density ratio m
 A = 1.0-(1.0-m)*eta^2;
 s = (r < eta);
 W = s.*(1.0-eta^2+m*(eta^2-r.^2))/A + (!s).*(1.0-r.^2)/A;
 dWdr = -2.0*(s*m+!s).*r/A;
end
