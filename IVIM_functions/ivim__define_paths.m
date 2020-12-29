function ivim__define_paths(iSubj,iPrePost)
%%
global vars init_vars


if iPrePost == 1
    fileID = fullfile(vars.dir.ResultsPath,[init_vars.Subj_list{iSubj},'_Pre_Results']);
elseif iPrePost == 2
    fileID = fullfile(vars.dir.ResultsPath,[init_vars.Subj_list{iSubj},'_Post_Results']);
elseif iPrePost == 3
    fileID = fullfile(vars.dir.ResultsPath,[init_vars.Subj_list{iSubj},'_Post_reg2Pre_Results']);
end

if exist([fileID,'.mat'])>0
    disp(['Results for ',init_vars.Subj_list{iSubj},' ',num2str(iPrePost),' already exists... '])
    vars.keepgoing = 0;
else
    vars.keepgoing = 1;
end

if vars.keepgoing == 1;
    if exist(vars.dir.ResultsPath)~=7
        mkdir(vars.dir.ResultsPath)
    end
    
    
%     if exist(fullfile(vars.dir.ResultsPath,'QA'))==7
% %         rmdir(fullfile(vars.dir.ResultsPath,'QA'),'s')
%     end
%     
    if exist(fullfile(vars.dir.QApath,'FitMean'))~=7
        mkdir(fullfile(vars.dir.QApath))
        mkdir(fullfile(vars.dir.QApath,'Maps'))
        mkdir(fullfile(vars.dir.QApath,'Maps','SliceMaps'))
        mkdir(fullfile(vars.dir.QApath,'FitMean'))
        mkdir(fullfile(vars.dir.QApath,'FitMean','BySlice'))
        mkdir(fullfile(vars.dir.QApath,'Presorting'))
        mkdir(fullfile(vars.dir.QApath,'T2'))
    end
end