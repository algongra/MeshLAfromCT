function indx = getlastzeroes(myarray)

indx=length(myarray);
while indx > 0
  if myarray(indx) == 0 
   indx = indx - 1;
  else
   break
  end
end
%indx = length(myarray) + indx;

return 
end

