function itemsByGroup = getcmplxgroups(self, ndim)


coords = {'x','y','z'};
imagre = {'im','re'};

% BUFFERS
buffers = {};
for j = 1:length(imagre)
   ff = imagre{j};
   myvar = ['p_' ff];
   buffers = [buffers myvar]; 
   for i = 1:3
      myvar = ['u' coords{i} '_' ff];
      buffers = [buffers myvar]; 
   end

end
%% MESHEU
%mesheu = {}; 
%for i = 1:ndim
%   myvar = [coords{i} 'p'];
%   mesheu = [mesheu myvar]; 
%   for j = 1:3
%      myvar = [coords{i} 'u' coords{j}];
%      mesheu = [mesheu myvar]; 
%   end
%   myvar = ['d' coords{i}];
%   mesheu = [mesheu myvar]; 
%end
%% MESHLAG
%meshlag = {}; 
%for i = 1:ndim
%   myvar = [coords{i} 'ib'];
%   myvar2 = ['f' coords{i}];
%   meshlag = [meshlag myvar myvar2]; 
%end
%
%% PARAMS
%params = {'Re','tf','geomfilename','nsdv','sdv',...
%          'testcase','resfilename','svnversion',...
%          'machine','makefile','ndim',...
%          'tsave','traw','tstats1samp','tstats1write', ...
%          'kz', 'bflowfilename'};
%for i = 1:ndim
%    myvar1 = ['n'  coords{i}];
%    myvar2 = ['nd' coords{i}];
%    myvar3 = [coords{i} '0'];
%    myvar4 = [coords{i} 'f'];
%    params = [params myvar1 myvar2 myvar3 myvar4];
%end

itemsByGroup.buffers = buffers;

return
end
