# Welcome to MeshLAfromCT

A library to automatically generate computational grids compatible with TUCAN from CT segmentations of left atria anatomies

## Installing MeshLAfromCT (libraries and software required to run it)

1. Generate binary files of iso2mesh (./iso2mesh/bin directory) following instructions in README.txt
2. Generate binary files of cpd2 (./ucsd-cvil_tools/cpd2/core/mex and ucsd-cvil_tools/cpd2/core/FGT directories) following instructions in INSTALL.txt
3. Check MATLAB is installed

## Using MeshLAfromCT
**1. Access** the **mfiles directory inside MeshLAfromCT directory** (this last one is the absolute path where you cloned this repository :p)
   ```
   cd ./install
   ```

**2. Copy seg2mesh_31174_10modes_nomres.m** (template file) **to seg2mesh_[NEW_CASE_ID]_[NUMBER_OF_FOURIER_MODES]modes_[SPATIAL_RESOUTION].m** file, **e.g.:**

   ```
   cp seg2mesh_31174_10modes_nomres.m seg2mesh_30701_10modes_nomres.m
   ```

**3. Modify inputs editing seg2mesh_[NEW_CASE_ID]_[NUMBER_OF_FOURIER_MODES]modes_[SPATIAL_RESOUTION].m** using your favourite editor, **e.g.:**

   ```
   vim seg2mesh_30701_10modes_nomres.m
   ```
   
   Example of modified inputs:
   ```
   .
   .
   .
   % Name for case, will be used to save/load files
   caseid = '30701'; % THIS INPUT WAS caseid = '31174'; in seg2mesh_31174_10modes_nomres.m

   % Extension of the segmentation files
   seg_ext = '.nii.gz';

   % times when segmentations were obtained
   seg_time = [17 27 40 54 67 80 94]/100; % THIS INPUT WAS seg_time = [0 7 14 21 28 36 43 50 57 64 71 78 85 92 99]/100; in seg2mesh_31174_10modes_nomres.m
   .
   .
   .
   % Inclusion tags
   tags.LA = 2;
   tags.LAA = 3;
   % Neighbouring tags
   tags.LV = 1;
   tags.PVs = 4:7; % THIS INPUT WAS tags.PVs = 5:7; in seg2mesh_31174_10modes_nomres.m
   
   % Cardiac period [s]
   T = 1;
   % Number of pulmunary veins
   nvein = 4; % THIS INPUT WAS nvein = 3; in seg2mesh_31174_10modes_nomres.m
   .
   .
   .
   mesh_keepratio = 0.4; % THIS INPUT WAS mesh_keepratio = 1; in seg2mesh_31174_10modes_nomres.m

   % Number of smoothing iterations on triangular mesh (with iso2mesh)
   smooth_iter = 2;

   % Max Area of faces [mm^2]
   MaxArea = 1.5*(10*dh).^2;

   % Multiplying factor of mean_area of mesh triangles as threshold
   % to detect open faces of PV and MV
   % [need to decrease with increasing resolution: 384 ~ 1, 256 ~ 0.5]
   % Default value: 0.5
   faces_cyl_thresh0 = 0.5; % THIS INPUT WAS faces_cyl_thresh0 = 0.95; in seg2mesh_31174_10modes_nomres.m 
   ```

**4. Run seg2mesh_[NEW_CASE_ID]_[NUMBER_OF_FOURIER_MODES]modes_[SPATIAL_RESOUTION].m** in MATLAB to generate computational grids compatible with TUCAN, **e.g.:**

   Running MATLAB in the background
   ```
   nohup matlab -nodisplay -r seg2mesh_30701_10modes_nomres > seg2mesh_30701_10modes_nomres.out 2> seg2mesh_30701_10modes_nomres.err &
   ```
   or open MATLAB:
   ```
   matlab -nodisplay
   ```
   and execute seg2mesh_30701_10modes_nomres.m in MATLAB Workspace:
   ```
   seg2mesh_30701_10modes_nomres
   ```
