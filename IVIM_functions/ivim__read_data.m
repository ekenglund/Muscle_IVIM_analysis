function ivim__read_data(iSubj,iPrePost)
global vars init_vars

% clear vars.ivim_info vars.mask_info vars.anat_info vars.AnatNames vars.mask_info_Mult vars.mask_info_ES
%%

% overwrite the init from previously
% if exist('/Volumes/MRI_Projects')
%     vars.dir.basedir = fullfile('/Volumes/MRI_Projects/3T/Spine/IVIM/Processed_Data'); % this is the base directory where the results, rois, and .nii files are located
%     vars.dir.base_dicomdir = fullfile('/Volumes/MRI/3T/Spine/IVIM/R03_Data'); % this is the base directory where the dicoms are located
% elseif exist('/Volumes/WardServer/MRI_Projects')
%     vars.dir.basedir = fullfile('/Volumes/WardServer/MRI_Projects/3T/Spine/IVIM/Processed_Data'); % this is the base directory where the results, rois, and .nii files are located
%     vars.dir.base_dicomdir = fullfile('/Volumes/WardServer/MRI/3T/Spine/IVIM/R03_Data'); % this is the base directory where the dicoms are located
% end
% vars.iPrePostMax = 3;
% vars.Subj_list = ...
%     {'HC01_01','HC01_02','HC03_01','HC03_02','HC04_01','HC04_02','HC05_01',...
%     'HC05_02','HC06_01','HC06_02','HC07_01','HC07_02','HC08_01','HC08_02',...
%     'HC09_01','HC10_01','HC11_01','HC11_02','HC12_01','HC12_02','HC13_01',...
%     'HC13_02','HC14_01','HC15_01','HC15_02','HC16_01','HC16_02','HC17_01'};
%%
if vars.MC_option == 0
    data_name_appendage = '.nii.gz';
elseif vars.MC_option == 1
    data_name_appendage = '.nii.gz';
elseif vars.MC_option == 2
    data_name_appendage = 'AFNIMC.nii.gz';
elseif vars.MC_option == 3
    data_name_appendage = 'noDC.nii';
%     data_name_appendage = 'manualDC.nii'; %
elseif vars.MC_option == 4
    data_name_appendage = 'manualDC_denoised.nii'; % change this back to 'manualDC_denoised.nii'
end

if iPrePost == 1
    disp(['Analyzing ',init_vars.Subj_list{iSubj},'_Pre ',data_name_appendage])
elseif iPrePost == 2
    disp(['Analyzing ',init_vars.Subj_list{iSubj},'_Post ',data_name_appendage])
elseif iPrePost == 3
    disp(['Analyzing ',init_vars.Subj_list{iSubj},'_PostReg2Pre ',data_name_appendage])
end


disp('Reading in data... ')

vars.dir.Subj_dir = dir(fullfile(vars.dir.basedir,[init_vars.Subj_list{iSubj},'*']));
vars.dir.Subj_dir = fullfile(vars.dir.Subj_dir.folder,vars.dir.Subj_dir.name);
vars.dir.dicomdir = dir(fullfile(vars.dir.base_dicomdir,[init_vars.Subj_list{iSubj},'*'],'e*'));
vars.dir.dicomdir = fullfile(vars.dir.dicomdir.folder,vars.dir.dicomdir.name);
vars.dir.datadir = dir(fullfile(vars.dir.Subj_dir,'Data'));
vars.dir.datadir = fullfile(vars.dir.datadir(1).folder);

roidir = fullfile(vars.dir.Subj_dir,'ROIs');
roidirs = dir(fullfile(roidir,'s*'));
roidirs_numbers(1) = str2num(roidirs(1).name(2:end));
roidirs_numbers(2) = str2num(roidirs(2).name(2:end));
resultsdir = fullfile(vars.dir.Subj_dir,'Results');
if iPrePost == 1
    dataname = 'IVIM_Pre*';
    vars.dir.MaskDir = fullfile(roidir,['s',num2str(min(roidirs_numbers))]);
    name_offset = 0;
elseif iPrePost == 2
    dataname = 'IVIM_Post*';
    vars.dir.MaskDir = fullfile(roidir,['s',num2str(max(roidirs_numbers))]);
    name_offset = 1;
elseif iPrePost == 3
    dataname = 'IVIM_Post*';
    data_name_appendage = '_reg2pre.nii*';
    vars.dir.MaskDir = fullfile(roidir,['s',num2str(min(roidirs_numbers))]);
    name_offset = 1;
end
fulldataname = dir(fullfile(vars.dir.datadir,[dataname,data_name_appendage]));
if isempty(fulldataname)
    disp(['No data found with ',data_name_appendage,'... Checking for noDC_AFNIMC_denoised.nii'])
    fulldataname = dir(fullfile(vars.dir.datadir,[dataname,'noDC_AFNIMC_denoised.nii']));
end
vars.fulldataname = fulldataname;

if length(vars.fulldataname)>1
    if vars.MC_option <2
        dataname_length = min(length(vars.fulldataname(1).name),length(vars.fulldataname(2).name));
        length_1 = length(fulldataname(1).name);
        length_2 = length(fulldataname(2).name);
        if dataname_length == length_1
            idx = 1;
        elseif dataname_length == length_2
            idx = 2;
        end
        vars.fulldataname = vars.fulldataname(idx).name;
    end
else
    vars.fulldataname = vars.fulldataname.name;
end
if vars.fulldataname(end)=='z'
    name = vars.fulldataname(1:end-7);
elseif vars.fulldataname(end)=='i'
    name = vars.fulldataname(1:end-4);
end

vars.dir.ResultsPath = fullfile(resultsdir,['20200214_',name,vars.dir.results_folder_name]);
vars.dir.QApath = fullfile(vars.dir.ResultsPath,'QA');

%%
if vars.fulldataname(end)=='z'
    raw_nii = load_untouch_nii(vars.fullfile(vars.dir.datadir,[name,'.nii.gz']));
elseif vars.fulldataname(end)=='i'
    raw_nii = load_untouch_nii(fullfile(vars.dir.datadir,[name,'.nii']));
end
b_images_raw = raw_nii.img;
b_images_raw = rot90(b_images_raw,-1);


if vars.MC_option == 3
    b_images_raw = reshape(b_images_raw,[128,128,22,150]);
end
    
    
vars.ivim.b_images_raw = b_images_raw(:,end:-1:1,:,:);


ImagePath_ES = fullfile(vars.dir.MaskDir,'ErectorSpinae');
ImagePath_Mult = fullfile(vars.dir.MaskDir,'Multifidus');

MaskNames_ES = dir(fullfile(ImagePath_ES,'*dcm*'));
MaskNames_Mult = dir(fullfile(ImagePath_Mult,'*dcm*'));

vars.dir.IVIMPath = fullfile(vars.dir.dicomdir,name(10+name_offset:15+name_offset));
a=dir(fullfile(vars.dir.IVIMPath,'*CFMRI*'));
if isfield(vars.ivim,'ivim_info')
    vars.ivim = rmfield(vars.ivim,'ivim_info');
end
for k = 1:vars.params.slices
    vars.ivim.ivim_info(k) = dicominfo(fullfile(vars.dir.IVIMPath,a(k).name));
end
clear a

ImagePathAnat = fullfile(vars.dir.dicomdir,vars.dir.MaskDir(end-5:end));
AnatNames = dir(fullfile(ImagePathAnat,'*CFMRI*'));

%% Read in masks and anatomical images
for k =1:vars.params.slices
    if size(dicomread(fullfile(ImagePath_Mult,MaskNames_Mult(k).name)),1) == 256
        mask_highres_Mult(:,:,k) = dicomread(fullfile(ImagePath_Mult,MaskNames_Mult(k).name));
        mask_info_Mult(k) = dicominfo(fullfile(ImagePath_Mult,MaskNames_Mult(k).name));
    else
        disp(['MASK IS INCORRECT FOR MULTIFIDUS SLICE NUMBER ',num2str(k)])
        mask_highres_Mult(:,:,k) = mask_highres(:,:,k-1);
        mask_info_Mult(k) = mask_info(k-1);
    end
end

for k =1:length(MaskNames_ES)
    if size(dicomread(fullfile(ImagePath_ES,MaskNames_ES(k).name)),1) == 256
        mask_highres_ES(:,:,k) = dicomread(fullfile(ImagePath_ES,MaskNames_ES(k).name));
        mask_info_ES(k) = dicominfo(fullfile(ImagePath_ES,MaskNames_ES(k).name));
    else
        disp(['MASK IS INCORRECT FOR SLICE NUMBER ',num2str(k)])
        mask_highres_ES(:,:,k) = vars.mask_highres(:,:,k-1);
        mask_info_ES(k) = mask_info(k-1);
    end
end
if isfield(vars,'mask')
    if isfield(vars.mask,'anat_info')
        vars.mask = rmfield(vars.mask,'anat_info');
    end
end

for k = 1:vars.params.slices
    vars.mask.anat_highres(:,:,k) = dicomread(fullfile(ImagePathAnat,AnatNames(k).name));
    vars.mask.anat_info(k) = dicominfo(fullfile(ImagePathAnat,AnatNames(k).name));
end

for k = 1:length(MaskNames_ES)
    vars.mask.mask_highres(:,:,k) = mask_highres_ES(:,:,k);
end
if length(MaskNames_ES)<vars.params.slices
    for k = length(MaskNames_ES)+1:vars.params.slices
        vars.mask.mask_highres(:,:,k) = mask_highres_Mult(:,:,k);
    end
end
vars.mask.MultMask_highres=mask_highres_Mult.*vars.mask.mask_highres;
vars.mask.ESMask_highres=vars.mask.mask_highres-vars.mask.MultMask_highres;



