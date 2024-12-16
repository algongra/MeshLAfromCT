function adjustfigures(figlist,varargin)
% @author M.Moriche
% @date 17-12-2013
%
% @brief function to adjust axis of several figures in order 
%        that they are comparable
% 
% @details
%
% Mandatory arguments:
% - figlist: list of figure handles to be adjusted.
% 
% Optional arguments:
% - [x|y][min|max]: [lower|upper] bound for [x|y] axis.
%
% Example:
%  @code
%  adjustfigures({1,2,3})
%  adjustfigures({1,2,3},'xmin',-1,'xmax',1)
%  adjustfigures({1,2,3},'ymin',0,'ymax',12.3)
%  @endcode
%

% default [x|y][min|max] calculated
xmax = -99999;
xmin = -xmax;
ymax = xmax;
ymin = xmin;
%
for i=1:length(figlist)
    fig = figlist{i};
    x = get(get(fig,'CurrentAxes'),'XLim');
    y = get(get(fig,'CurrentAxes'),'YLim');
    if x(1) < xmin
        xmin = x(1)
    end
    if x(2) > xmax
        xmax = x(2)
    end
    if y(1) < ymin
        ymin = y(1)
    end
    if y(2) > ymax
        ymax = y(2)
    end
end

misc.assigndefaults(varargin{:});
%
for i=1:length(figlist)
    fig = figlist{i};
    set(get(fig,'CurrentAxes'),'XLim',[xmin xmax]);
    set(get(fig,'CurrentAxes'),'YLim',[ymin ymax]);
end
%
%
return
end
