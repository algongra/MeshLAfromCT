% A wrapper to cpd_register, allowing the registration of wavefront obj
% files.
%
% Parameters:
% - f_fn: Fixed file name.
% - m_fn: Moving file name.
% - o_fn: Output file name.
% - opt: CPD options.
%
% Relies upon:
% - read_wobj
% - write_wobj
% - cpd_register

function cpd_register_obj(f_fn, m_fn, o_fn, opt)
    % Add repo/tools to your path.
    addpath(genpath('C:\Users\woodford\ucsd_cvi_repo\tools'))
    addpath(genpath('C:\Users\woodford\ucsd_cvi_repo\squeez\tools\cpd2\core\'))

    % Read examples.
    f_obj = read_wobj(f_fn);
    m_obj = read_wobj(m_fn);
    o_obj = m_obj;

    % Perform the CPD registration.
    [Transform, ~]=cpd_register(f_obj.vertices, m_obj.vertices, opt);

    o_obj.vertices = Transform.Y;

    write_wobj(o_obj, o_fn)
end