function data = getplotinfo(fig)
% @author M.Moriche
% @date 17-12-2013 by M.Moriche \n
%       Documented
%
% @brief Function to extract information about a figure
%
% @details
%
% The information looked for is how data is represented
%
% Output
% - data: cell of structs where every struct represent a series
%         of the plot.
% 
% Mandatory arguments
% - fig: figure handle
%  
% No optional arguments
%

series = get(fig,'Children');
series = series(1);
series = get(series,'Children');

data = {};
for iser = 1:length(series)
   se = get(series(iser));

   sedata.LineStyle   = se.LineStyle;
   sedata.Color       = se.Color    ;
   sedata.DisplayName = se.DisplayName;
   sedata.Marker      = se.Marker;

   data = [data sedata];
end

end
