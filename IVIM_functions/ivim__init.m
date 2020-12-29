function ivim__init
global vars init_vars

% vars.fitMaps = 1;
% vars.fitBySlices = 1; 
% vars.fitByTime = 1;
% vars.fitByTimeAndSlice = 1;

vars.params.tensor = 5;
vars.params.dirs = 149;
vars.params.highb_threshold = 200;
vars.params.read = 128;     % Number of dummy_scans before OxFlow sequence (8 dummy scans)
vars.params.pe = 128;              % Number of lines in phase encoding direction (208 lines)
vars.params.slices = 22;         % Number of images acquired (4 interleaves)

if vars.saveQA ~=1
    input('Are you sure you do not want to save QA figures?')
end
if vars.vis == 1
    input('Are you sure you want to keep all QA figures open?')
end
% if vars.fitMaps ~= 1 
%     input('Are you sure you do not want to fit IVIM parameter maps?')
% end
% if vars.fitBySlices ~= 1 
%     input('Are you sure you do not want to fit IVIM parameter by slice?')
% end
% if vars.fitByTime ~= 1 
%     input('Are you sure you do not want to fit IVIM parameter by time?')
% end
% if vars.fitByTimeAndSlice ~= 1 
%     input('Are you sure you do not want to fit IVIM parameter by time and slice?')
% end
if exist('/Volumes/MRI_Projects')
    vars.dir.basedir = fullfile('/Volumes/MRI_Projects/3T/Spine/IVIM/Processed_Data'); % this is the base directory where the results, rois, and .nii files are located
    vars.dir.base_dicomdir = fullfile('/Volumes/MRI/3T/Spine/IVIM/R03_Data'); % this is the base directory where the dicoms are located
elseif exist('/Volumes/WardServer/MRI_Projects')
    vars.dir.basedir = fullfile('/Volumes/WardServer/MRI_Projects/3T/Spine/IVIM/Processed_Data'); % this is the base directory where the results, rois, and .nii files are located
    vars.dir.base_dicomdir = fullfile('/Volumes/WardServer/MRI/3T/Spine/IVIM/R03_Data'); % this is the base directory where the dicoms are located
end

init_vars.Subj_list = ...
    {'HC01_01','HC03_02','HC04_01','HC05_01','HC06_01','HC07_01','HC08_01',...
    'HC09_01','HC10_01','HC11_01','HC12_01','HC13_01','HC14_01','HC15_01',...
    'HC16_01','HC17_01','LS01_01','LS02_01','LS03_01','LS04_01','LS05_01',...
    'LS06_01','LS07_01','LS08_01','LS09_01','LS10_01','LS11_01','LS12_01',...
    'LS13_01','LS14_01','LS15_01','LS16_01','LS17_01','LS18_01',... 
    'HC01_02','HC03_01','HC04_02','HC05_02','HC06_02','HC07_02','HC08_02',...
    'HC11_02','HC12_02','HC13_02','HC15_02','HC16_02',...
    'LS03_02','LS04_02','LS10_02','LS14_02','LS14_03'};

 
%     
vars.MC_option = 4; % 0 = no MC, 1 = 0 = no motion correction, 1 = matlab motion correction, 2 = afni motion correction, 3 = + manual DC, 4 = + manaul DC + denoised

if vars.MC_option == 0
    vars.dir.results_folder_name = '_Results_noMC';
    vars.iPrePostMax = 2;
elseif vars.MC_option == 1
    vars.dir.results_folder_name = '_Results';
    vars.iPrePostMax = 2;
elseif vars.MC_option == 2
    vars.dir.results_folder_name = '_Results';
    vars.iPrePostMax = 3;
elseif vars.MC_option == 3
    vars.dir.results_folder_name = '_Results';
    vars.iPrePostMax = 2;
elseif vars.MC_option == 4
    vars.dir.results_folder_name = '_Results';
    vars.iPrePostMax = 3;
end



