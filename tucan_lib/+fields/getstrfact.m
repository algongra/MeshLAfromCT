function k = getstrfact(dom,domr,dr,kmin,kmax)

% @author A.Gonzalo
% 
% @brief function to find out the stretching factor needed to get a mesh
%        with boundary conditions fixed at a spatial domain
% 
% @date  03-03-2015 by A.Gonzalo
%        Created
%
% @details
%
%  MANDATORY ARGUMENTS:
%  --------------------
%
%  - dom: initial/final point of physical domain
%
%  - domr: intial/final physical point of refined section
%
%  - dr: distance between points in refined section
%
%  - kmin: minimum stretching factor allowed
%
%  - kmax: maximum stretching factor allowed
%
%  OUTPUT ARGUMENTS:
%  -----------------
%
%  - k: stretching factor
%
% @code
% k = fields.getstrfact(dom,domr,dr,kmin,kmax);
% @endcode

% get maximum stretching factor
sfmax = 1 + kmax;
% get minimum stretching factor
sfmin = 1 + kmin;
% get stretching section distance
S = abs(dom-domr);
% get maximum number of points
n = floor(S/dr);

% n must be odd to get at least one real root 
%if mod(n,2) == 0
%   n = n-1;
%end

% coefficients of the polynomial
C = S/dr;
B = (1+C);
A = 1;

% initialize k
sf = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% needed in displays
iter = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while sf<=sfmin || sf>=sfmax
      % create the polynomial
      clear p; p = zeros(1,n+2);
      p(1) = A; p(end-1) = -B; p(end) = C;
      
      % get roots of the polynomial
      myRoots = roots(p);

      % get real roots of the polynomial
      realindx = myRoots == real(myRoots);
      myRealRoots = myRoots(realindx);
      
      % We are looking for a a stretching factor greater than 1
      tol = 10e-12;
      posindx = myRealRoots > (1+tol);
      
      sf = myRealRoots(posindx);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % displays
      iter = iter + 1;
      fprintf('\niteration %i\n',iter)
      fprintf('n = %i\n',n)
      fprintf('sf = %17.15f\n',sf)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if sf > sfmax
         n = floor(n+1);
         %n = floor(3*n/2);
      elseif sf < sfmin
         n = floor(n-1);
         %n = floor(n/2);
      end
      
      % get possible stretching factors
      cond1 = myRealRoots > sfmin; cond2 = myRealRoots < sfmax;
      sfindx = logical(cond1.*cond2);
      
      sf = myRealRoots(sfindx);
      
      % n must be odd to get at least one real root 
      %if mod(n,2) == 0
      %   n = n-1;
      %end
      
      % checks
      if isempty(sf) == 1
         sf = 0;
      end
      if n<3
         k = [];
         aggfrmt = ['It is not possible to find a stretching factor '...
                    'between kmin and kmax with those values of dom, '...
                    'domr and dr'];
         fprintf('\n%s\n\n',aggfrmt)
         return
      end
end

k = sf - 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% more displays
u(1) = domr;
i = 1;
j = 1;
while n >= i
      i = i + 1;
      u(i) = u(i-1) + dr*min(sf)^j;
      j = j + 1;
end
fprintf('\nFirst velocity point would be located at %17.15f\n',u(1))
fprintf('Last velocity point would be located at %17.15f\n\n',u(end))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
