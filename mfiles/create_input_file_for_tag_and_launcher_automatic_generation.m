function create_input_file_for_tag_and_launcher_automatic_generation(pth,keys,values);

% Assign values to keys
for i = 1:length(values)
    value = values{i};
    
    % Check if the value is a number or string
    if isnumeric(value)
        % If it's a number, assign it directly
        eval([keys{i} ' = ' num2str(value) ';']);
    elseif ischar(value) || isstring(value)
        % If it's a string, handle it accordingly
        eval([keys{i} ' = ''' value ''';']);
    end
end

% Absolute path of the file
basenm = sprintf('tucanHGPU_preprocessing_input_case%s_%imodes_%s.sh',...
                 caseid,nmodes,resnm);
flnm = fullfile(pth,basenm);

% Write file
%
% Open file
fid = fopen(flnm,'w');
%
% Check if the file opened successfully
if fid == -1
   error('Failed to open the file for writing');
end
%
% Write file content line by line
fprintf(fid,'#!/bin/bash\n\n');
fprintf(fid,'################################################################################\n');
fprintf(fid,'#                                    INPUTS                                    #\n');
fprintf(fid,'################################################################################\n\n');
fprintf(fid,'# Paths definition:\n');
fprintf(fid,'# -----------------\n');
fprintf(fid,'# IBcode GPU repository absolute path\n');
fprintf(fid,'IBcodeGPU_repos=\n');
fprintf(fid,'#\n');
fprintf(fid,'# miniconda3 bin directory absolute path\n');
fprintf(fid,'minicondabindir=\n');
fprintf(fid,'#\n');
fprintf(fid,'# Mesher absolute path\n');
fprintf(fid,'mesher_dir=\n');
fprintf(fid,'#\n');
fprintf(fid,'# Directory where tucanHGPU simulation input files will be stored\n');
fprintf(fid,'input_files_dir=\n');
fprintf(fid,'#\n');
fprintf(fid,'# Directory where tucanHGPU tag to run this case_ID simulation will be created\n');
fprintf(fid,'tag_dir=\n');
fprintf(fid,'#\n');
fprintf(fid,'# Directory where tucanHGPU launcher to run this case_ID simulation will be created\n');
fprintf(fid,'launcher_dir=\n');
fprintf(fid,'#\n');
fprintf(fid,'# Directory where tucanHGPU simulation outputs will be stored\n');
fprintf(fid,'out_dir=\n');
fprintf(fid,'#\n');
fprintf(fid,'# Name definition:\n');
fprintf(fid,'# -----------------\n');
fprintf(fid,'# Name of IBcode GPU branch to be used\n');
fprintf(fid,'branch_name=tucanHNNgpu_RT2dw\n');
fprintf(fid,'#\n');
fprintf(fid,'# Case ID\n');
fprintf(fid,'case_ID=%s\n',caseid);
fprintf(fid,'#\n');
fprintf(fid,'# Number of Fourier modes used to describe LA endocardium motion\n');
fprintf(fid,'nmodes=%i\n',nmodes);
fprintf(fid,'#\n');
fprintf(fid,'# Eulerian spatial resolution name\n');
fprintf(fid,'# lowres --> low resolution (144x144x144)\n');
fprintf(fid,'# nomres --> nominal resolution (256x256x256)\n');
fprintf(fid,'resolution=%s\n',resnm);
fprintf(fid,'#\n');
fprintf(fid,'# Tag name (new tag created for this simulation run)\n');
fprintf(fid,'tag_name=${branch_name}_case${case_ID}_${nmodes}modes_${resolution}\n');
fprintf(fid,'#\n');
fprintf(fid,'# CFD simulation parameters:\n');
fprintf(fid,'# --------------------------\n');
fprintf(fid,'#  - Memory managing\n');
fprintf(fid,'# Buffers precision (''single'' or ''double''). Default value: ''double''\n');
fprintf(fid,'pre=''double''\n');
fprintf(fid,'#  - Parallelization\n');
fprintf(fid,'# Threads per block (Eulerian mesh). Default value: 6\n');
fprintf(fid,'TPB=6\n');
fprintf(fid,'# Threads per block (Lagrangian mesh). Default value: 181\n');
fprintf(fid,'TPBlag=181\n');
fprintf(fid,'#  - Simulation parameters\n');
fprintf(fid,'# Reynolds number. Default value: 25\n');
fprintf(fid,'Re=25\n');
fprintf(fid,'# Volumetric forcing in x, y and z directions. Default values: (0 0 0)\n');
fprintf(fid,'bodyf=(0 0 0)\n');
fprintf(fid,'# Number of veins. Default value: 4\n');
fprintf(fid,'nvein=%i\n',nvein);
fprintf(fid,'# Index body to identify first cork in Lagrangian mesh. Default value: 2+2*nvein\n');
fprintf(fid,'icork=$((2+2*${nvein}))\n');
fprintf(fid,'# Forcing relaxation used to impose BCs of velocity on corks. Default value: 0.1\n');
fprintf(fid,'bet=0.1\n');
fprintf(fid,'# Forcing relaxation used to impose BCs of RT on corks. Default value: 0.1\n');
fprintf(fid,'alp=0.1\n');
fprintf(fid,'# Forcing relaxation used to impose BCs of RT outside LA. Default value: 0.01\n');
fprintf(fid,'dtforc=0.01\n');
fprintf(fid,'# Hematrocrit value.\n');
fprintf(fid,'# NO DEFAULT (INPUT REQUIRED)\n');
fprintf(fid,'Hct=\n');
fprintf(fid,'# Non-Newtonian rheology model kinematic viscosity parameters [cm^2/s]\n');
fprintf(fid,'# NO DEFAULT (INPUTS REQUIRED)\n');
fprintf(fid,'nuzero=0.64\n');
fprintf(fid,'nuinf=0.04\n');
fprintf(fid,'# Non-Newtonian rheology activation based on RT (thixotropic blood behavior).\n');
fprintf(fid,'# Default values RTmean: 3, RTwidth: 2\n');
fprintf(fid,'RTmean=3\n');
fprintf(fid,'RTwidth=2\n');
fprintf(fid,'# Eulerian mesh resolution per box length. Default value: 256\n');
fprintf(fid,'ppL=%i\n',res);
fprintf(fid,'# Eulerian mesh box length in cm. Default value: 13\n');
fprintf(fid,'Lc=%g\n',Lc);
fprintf(fid,'# Eulerian mesh physical boundaies definition in cm.\n');
fprintf(fid,'# Default values: (-6.5 6.5 -6.5 6.5 -6.5 6.5)\n');
fprintf(fid,'# [i.e, xb = -6.5, xf = 6.5, yb = -6.5, yf = 6.5, zb = -6.5, zf = 6.5]\n');
fprintf(fid,'xyzbox=(%g %g %g %g %g %g)\n',x0,xf,y0,yf,z0,zf);
fprintf(fid,'# dt used in CFD simulation. Default value: 5e-5\n');
fprintf(fid,'dt=5e-5\n');
fprintf(fid,'# Number of time steps run in the CFD simulation.\n');
fprintf(fid,'# Default value: 20000*20 # 20000 steps/cycle * 20 cycles\n');
fprintf(fid,'nsteps=\n');
fprintf(fid,'# Number of steps between frame saving. Default value: 500\n');
fprintf(fid,'nsave=\n');
fprintf(fid,'# Directory where tucanHGPU simulation outputs will be stored.\n');
fprintf(fid,'# NO DEFAULT (INPUT REQUIRED)\n');
fprintf(fid,'initpath=${out_dir}\n');
fprintf(fid,'# Name of simulation/case to be run.\n');
fprintf(fid,'# NO DEFAULT (INPUT REQUIRED)\n');
fprintf(fid,'basename=case${case_ID}_dw_${nmodes}modes_${resolution}_RTnonnewtonian1\n');
fprintf(fid,'# Restart flag option\n');
fprintf(fid,'# NO DEFAULT (INPUT REQUIRED)\n');
fprintf(fid,'# 0 --> All buffers (velocity and RT) start from initial conditions (all 0s)\n');
fprintf(fid,'#       iframe and istep0 start from 0\n');
fprintf(fid,'# 1 --> Velocities restart from file, RT starts initial conditions (all 0s)\n');
fprintf(fid,'#       iframe and istep0 start from 0\n');
fprintf(fid,'# 2 --> All buffers restart from file\n');
fprintf(fid,'#       iframe and istep0 start are read from resfilename\n');
fprintf(fid,'# 3 --> Velocities restart from file, RT starts initial conditions (all 0s)\n');
fprintf(fid,'#       iframe and istep0 start are read from resfilename\n');
fprintf(fid,'tocontinue=0\n');
fprintf(fid,'# Name of restart file\n');
fprintf(fid,'# NO DEFAULT (INPUT REQUIRED IF tocontinue NOT EQUAL TO ZERO)\n');
fprintf(fid,'resfilename=\n');
fprintf(fid,'#\n');
fprintf(fid,'# Ask before changing any of the following parameters:\n');
fprintf(fid,'#\n');
fprintf(fid,'# Frequence of cardiac cycle\n');
fprintf(fid,'freq=1\n');
fprintf(fid,'# Time sampling used to store additional inputs information\n');
fprintf(fid,'dtsamp=%g\n',dt_spln);
fprintf(fid,'# dt used to create pre-processing files\n');
fprintf(fid,'dtpre=5e-3\n');
fprintf(fid,'#\n');
fprintf(fid,'# Launcher header variables definition:\n');
fprintf(fid,'# -------------------------------------\n');
fprintf(fid,'# Klone group name\n');
fprintf(fid,'groupname=mambrino\n');
fprintf(fid,'# Klone resource partition to be used\n');
fprintf(fid,'partition=ckpt\n');
fprintf(fid,'# Simulation wall time format: DD-HH:MM:SS or HH:MM:SS)\n');
fprintf(fid,'runwalltime=30:00:00\n');
fprintf(fid,'# Total RAM (CPU) memory to be allocated for the job (format: 80G or 100M)\n');
fprintf(fid,'RAMmem=80G\n');
fprintf(fid,'# Generic resources [e.g. GPUs] (format: gpu:a100:1 or gpu:a40:1)\n');
fprintf(fid,'gputypeandnum=gpu:a100:1\n');
fprintf(fid,'# Number of CPUs allocated per task\n');
fprintf(fid,'nprocpertask=7\n');
fprintf(fid,'# Restart flag after simulation start running (with flag tocontinue\n');
fprintf(fid,'# 1 --> Velocities restart from file, RT starts initial conditions (all 0s)\n');
fprintf(fid,'#       iframe and istep0 start from 0\n');
fprintf(fid,'# 2 --> All buffers restart from file\n');
fprintf(fid,'#       iframe and istep0 start are read from resfilename\n');
fprintf(fid,'# 3 --> Velocities restart from file, RT starts initial conditions (all 0s)\n');
fprintf(fid,'#       iframe and istep0 start are read from resfilename\n');
fprintf(fid,'tocontinuelauncher=2\n');
fprintf(fid,'# Name of conda environment to be loaded\n');
fprintf(fid,'envname=gpu\n');
fprintf(fid,'# Launcher file template absolute path (including .sh extension)\n');
fprintf(fid,'template_launcher_file=${IBcodeGPU_repos}/lib/bashfiles/launcher_file_template.sh\n\n');
fprintf(fid,'################################################################################\n');
%
% Close the file after writing
fclose(fid);

% Change permission of exceution for user, groups, and all
system(sprintf('chmod 755 %s',flnm));

% Print message if it success
disp(sprintf('\nInputs files for tag and launcher generation saved in:\n %s',...
             flnm));

return

end
