function lns = text2lines(text,w)

lns = {};

i1 = 1;
i2 = w;
while length(text) > w
   i2 = min(w,length(text));
   lns = [lns {text(i1:i2)}];
   text = text(i2+1:end);
end
i2 = min(w,length(text));
lns = [lns {text(i1:i2)}];


return
end

