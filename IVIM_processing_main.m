%% Streamline IVIM Processing for a variety of different diffusion directions, number of repeats, etc
% 
close all,clear all,clc
% Depending on whether you're processing on the iMac or on your laptop, the paths will differ
% This next step sets the appropriate paths so you don't have to worry about it
if exist('/Volumes/Code') 
    addpath /Volumes/Code/m-files/
    addpath /Volumes/Code/m-files/NIfTI/
    addpath /Volumes/Code/m-files/afni_matlab/matlab/
    addpath /Volumes/Code/m-files/20130227_xlwrite/20130227_xlwrite/
    addpath /Volumes/Code/IVIM/IVIM_functions/
elseif exist('/Volumes/WardServer/Code')
    addpath /Volumes/WardServer/Code/m-files/
    addpath /Volumes/WardServer/Code/m-files/NIfTI/
    addpath /Volumes/WardServer/Code/m-files/afni_matlab/matlab/
    addpath /Volumes/WardServer/Code/m-files/20130227_xlwrite/20130227_xlwrite/
    addpath /Volumes/WardServer/Code/IVIM/IVIM_functions/
end
tic

global vars init_vars
%%

vars.vis = 0; % vis = 1 keeps all figures open
vars.saveQA = 1; % saveQA = 1 saves QA figures
vars.saveQAall = 0; % saveQA = 1 saves QA figures

% initialize IVIM processing paths and function
ivim__init % here is where you set your subject list, paths, etc.
%
for iSubj = 1:length(init_vars.Subj_list)% go through each subject in Subj_list
    for iPrePost = 1:vars.iPrePostMax
        close all
        ivim__read_data(iSubj,iPrePost)
        ivim__define_paths(iSubj,iPrePost)
        if vars.keepgoing == 1
            ivim__preprocessing % Median filter (figures 100s)
            ivim__mask_manipulation(iSubj,iPrePost) % (figures 110s)
            ivim__b_images_manipulation % (figures 120s)
            ivim__relativeT2 % (figures 130s)
            ivim__separate_Sbover1(iSubj,iPrePost)
            ivim__fit_map(iSubj,iPrePost) % (figures 160s)
            ivim__save_and_export_rmalloutliers_2(iSubj,iPrePost)
            ivim__fit_resolve_time_varymask(iSubj,iPrePost) % (figures 140s)
            ivim__fit_resolve_slice_varymask(iSubj,iPrePost) % (figures 150s)
            ivim__save_and_export(iSubj,iPrePost)
            if iPrePost == 1
                save(fullfile(vars.dir.ResultsPath,[vars.Subj_list{iSubj},'_Pre_Results']))
            elseif iPrePost == 2
                save(fullfile(vars.dir.ResultsPath,[vars.Subj_list{iSubj},'_Post_Results']))
            elseif iPrePost == 3
                save(fullfile(vars.dir.ResultsPath,[vars.Subj_list{iSubj},'_Post_reg2Pre_Results']))
            end
        end
        ivim__briksave(iSubj,iPrePost)
    end
end
toc
