function warningmessage(mystring)


display(' ')
display(' ')
myhigh = '*';
for i0 = 1:(80-1)
   myhigh = [myhigh '*'];
end
display(sprintf('%s\n',myhigh));
display(' ')
display(' ')
mylist = misc.text2lines(mystring, 80);
for i0 = 1:length(mylist)
   mys = mylist{i0};
   display(sprintf('%s\n', mys));
end
display(' ')
display(' ')
display(' ')
display(sprintf('%s\n',myhigh));
display(' ')

return
end
