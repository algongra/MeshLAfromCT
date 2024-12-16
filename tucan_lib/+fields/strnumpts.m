function [ns leng] = strnumpts(lstop,dstop,dr,sr)

% set a tolerance
tol = 10e-12;

% Get ncs(iic,iipos) and leng
% Get reference left/down/bottom|right/up/upper stretching
%  zone size and and length of those stretching zone
leng = 0;
lengold = 0;
delta = dr;
ns = 0;
nsold = 0;
while leng<=(lstop+tol) && delta<(dstop+tol)
       nsold = ns;
       ns = ns + 1;
       %
       lengold = leng;
       leng = leng + dr*sr^ns;
       %
%       deltaold = delta;
       delta = delta*sr;
end
ns = nsold;
leng = lengold;

end
