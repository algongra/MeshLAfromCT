function fig = plotlegend(data)
% @author M.Moriche
% @date 17-12-2013 by M.Moriche \n
%       Documented
%
% @brief Function to plot the legend item from the struct
%        data with information about how the series is 
%        represented
%
% @details
%
% The name of the series is the one saved in the field 
% DisplayName, this has been kept from Matlab rules.
% 
% If DisplayName is empty, therefore it does not need a 
% legend item, returns -1.
%
% Output
% - fig: figure handle of the legend plot.
%        Returns -1 if the data has no name to diplay
% 
% Mandatory arguments
% - data: struct with data of the series representation.
%
% No optional arguments
%

x0=0;
xf=1;
y0=0;
yf=0;

mysize = [0.3 0.1];
vv = [0 0 mysize];

conf.Units = 'inches';
conf.Position = [0 0 mysize];
conf.PaperPosition = [0 0 mysize];
conf.PaperSize     = mysize;
conf.Visible = 'off';

axconf.Units='normalized';
axconf.Position=[0.1 0.1 0.8 0.8];
axconf.Visible = 'off';

if ~isempty(data.DisplayName)
   fig = figure(conf);
   plot([x0 xf], [y0 yf], data);
   set(gca,axconf);
else
   fig = -1;
end

end
