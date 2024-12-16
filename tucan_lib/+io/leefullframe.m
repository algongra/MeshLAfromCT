function mystruct = leefullframe(basename,iext,varargin)
%fr = leefullframe(basename,iext,varlooked)
%fr = leefullframe(basename,iext,varlooked,path)
%
%  example fr = leefullframe(basename,3);
%  example fr = leefullframe(basename,3,{'ux','xux','yux'});
%  example fr = leefullframe(basename,3,{'ux','xux','yux'},'/data/');
items='*';
path=getenv('IBCODEDATA');
misc.assigndefaults(varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ext = sprintf('%04d.000',iext);
%[mypath '/' basename '/' basename '.' ext]
fnm = [path '/' basename '/' basename '.' ext];
mystruct = io.leefulldata(fnm,items);

return 
end
