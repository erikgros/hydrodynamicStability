function [I] = integ( fk, r )
% integral of a function, defined by its values fk at the Chebyshev nodes, over r
 N = length( fk );
 if (!isfinite(fk(end)))
  fk = fk(1:end-1);
  xk = sin(pi*[N-1:-2:1-N+2]'/(2*(N-1)));    % Chebyshev points on (-1, 1]
 else
  xk = sin(pi*[N-1:-2:1-N]'/(2*(N-1)));    % Chebyshev points on [-1, 1]
 end
 % using vander to find interpolating polynomial:
 A = vander( xk );
 pol = A \ fk;

 Ipol = polyint( pol ); % (polynomial integration)
 a = min( r );
 b = max( r );
 I = diff(polyval(Ipol,[-1 1]));
 I *= (b-a) * 0.5;
end
