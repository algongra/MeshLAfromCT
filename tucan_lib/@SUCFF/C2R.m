function C2R_z(self, beta, z, varargin);
%
% @author M.Moriche
% @brief Function to generate real flow field in z direction
%
% @details
%
% -MANDATORY ARGUMENTS: 
%   beta: wavenumber
%   z: z position to get the real values
%
% -OPTIONAL ARGUMENTS: 
%   varlist = {'ux', 'uy', 'uz', 'p'}
%
%

varlist = {'ux','uy','uz','p'};
misc.assigndefaults(varargin{:});

for i0 = 1:length(varlist)
   varnm = varlist{i0};
   %var = real(self.([varnm '_re'])*cos(beta*z) + j*self.([varnm '_im'])*sin(beta*z));
   var = self.([varnm '_re'])*cos(beta*z) - self.([varnm '_im'])*sin(beta*z);
   self.setfield(varnm, var);
end

return
end
