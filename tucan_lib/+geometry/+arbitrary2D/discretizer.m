function xy = discretizer( mygeomfun, s0, sf, varargin)

maxds = (sf-s0)/10;
maxangle = 0.005;
n0 = 10000;
misc.assigndefaults(varargin{:});


xy = zeros(n0,2);

ds = (sf-s0)/100;

s = s0;
xynow = mygeomfun(s);

i = 1;
xy(1,:) = xynow;
i = i + 1;

while s < sf

   flag = false;
   n2 = 0;
   ds = ds*2;

   while ~flag

      ds = min([ds/2, sf-s, maxds]);

      % Trial 1
      s1 = s + ds;
      xy1 = mygeomfun(s1);

      % Trial 2
      s2 = s + 0.1*ds;
      xy2 = mygeomfun(s2);

      % check step
      flag = geometry.arbitrary2D.checkstep(xynow, xy1, xy2, maxds, maxangle);

      n2 = n2 + 1;

   end

   xynow = xy1;

   xy(i,:) = xynow;
   i = i + 1;
   s = s + ds;

   % if it went smooth, increase ds
   if n2 == 1
      ds = ds*2;
   end

end

xy = xy(1:(i-1),:);

return
end
