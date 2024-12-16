function [facecell varargout]=finddisconnsurf2(f)
%
% facecell=finddisconnsurf(f)
% 
% subroutine to extract disconnected surfaces from a 
% cluster of surfaces
% 
% author: Qianqian Fang (fangq@nmr.mgh.harvard.edu)
% Date: 2008/03/06
% author: Alejandro Gonzalo (algongra@uw.edu)
% Date: 2023/05/15 --> Modified to add an optional output
%
% input: 
%     f: faces defined by node indices for all surface triangles
%
% output:
%     facecell: separated disconnected surface node indices
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%
f0 = f;

faceid = {};
facecell={};
subset=[];
while(~isempty(f))
	idx0=reshape(ismember(f0,f0(1,:)), size(f0));
	idx=reshape(ismember(f,f(1,:)), size(f));
	ii=find(sum(idx,2));
	while(~isempty(ii))
		if(isempty(ii)) break; end
		%ii=unique(ii);
		subset(end+1:end+length(ii),:)=f(ii,:);
		f(ii,:)=[];
		idx0=reshape(ismember(f0,subset), size(f0));
		idx=reshape(ismember(f,subset), size(f));
	   ii=find(sum(idx,2));
	end
	if(~isempty(subset))
      ii0=find(sum(idx0,2));
      faceid{end+1} = ii0;
		facecell{end+1}=subset;
      subset=[];
	end
end

if nargout == 2
   varargout{1} = faceid;
end
