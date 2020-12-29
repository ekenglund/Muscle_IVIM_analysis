function ivim__read_data(iSubj,iPrePost)
global vars
%%
disp('Saving BRIK files...')
% 


if exist('/Volumes/MRI_Projects') && strcmp(vars.dir.datadir(1:20),'/Volumes/MRI_Project')
    datadir = vars.dir.datadir;
    ResultsPath = vars.dir.ResultsPath;
elseif exist('/Volumes/MRI_Projects') && strcmp(vars.dir.datadir(1:20),'/Volumes/WardServer/')
    datadir = ['/Volumes/',vars.dir.datadir(21:end)];
    ResultsPath = ['/Volumes/',vars.dir.ResultsPath(21:end)];
elseif exist('/Volumes/WardServer/MRI_Projects') && strcmp(vars.dir.datadir(1:20),'/Volumes/WardServer/')
    datadir = vars.dir.datadir;
    ResultsPath = vars.dir.ResultsPath;
elseif exist('/Volumes/WardServer/MRI_Projects') && strcmp(vars.dir.datadir(1:20),'/Volumes/MRI_Project')
    datadir = ['/Volumes/WardServer/',vars.dir.datadir(10:end)];
    ResultsPath = ['/Volumes/WardServer/',vars.dir.ResultsPath(10:end)];
end


if ~exist(fullfile(ResultsPath,'QA','AFNI'))
    mkdir(fullfile(ResultsPath,'QA','AFNI'))
    AFNI_exist = 0;
else 
    AFNI_exist = 1;
    disp('BRIK file already exist...')
end
cd(fullfile(ResultsPath,'QA','AFNI'))


if vars.keepgoing == 0 && AFNI_exist == 0;
    a = dir(fullfile(ResultsPath,'*_Results.mat'));
    load(fullfile(ResultsPath,a.name));
end
    
if iPrePost == 1
    if exist(fullfile(datadir,'PreT1+orig.BRIK'))==2
        T1_exist = 1;
        T1 = BrikLoad(fullfile(datadir,'PreT1+orig.BRIK'));
        [err,T1_header] = BrikInfo(fullfile(datadir,'PreT1+orig.HEAD'));
    else
        display(['Need PreT1 for ',vars.Subj_list{iSubj}])
        T1_exist = 0;
    end
elseif iPrePost == 2
    if exist(fullfile(datadir,'PostT1+orig.BRIK'))==2
        T1_exist = 1;
        T1 = BrikLoad(fullfile(datadir,'PostT1+orig.BRIK'));
        [err,T1_header] = BrikInfo(fullfile(datadir,'PostT1+orig.HEAD'));
    else
        display(['Need PostT1 for ',vars.Subj_list{iSubj}])
        T1_exist = 0;
    end
elseif iPrePost == 3
    if exist(fullfile(datadir,'PostT1reg+orig.BRIK'))==2
        T1_exist = 1;
        T1 = BrikLoad(fullfile(datadir,'PostT1reg+orig.BRIK'));
        [err,T1_header] = BrikInfo(fullfile(datadir,'PostT1reg+orig.HEAD'));
    else
        display(['Need PostT1reg for ',vars.Subj_list{iSubj}])
        T1_exist = 0;
    end
end

% 
% if exist(fullfile(ResultsPath,'QA','AFNI'))~=7
%     mkdir(fullfile(ResultsPath,'QA','AFNI'))
% end
% 

if ~exist(fullfile(ResultsPath,'QA','AFNI'))
    mkdir(fullfile(ResultsPath,'QA','AFNI'))
    AFNI_exist = 0;
else 
    AFNI_exist = 1;
    disp('BRIK file already exist...')
end
cd(fullfile(ResultsPath,'QA','AFNI'))
    
if T1_exist == 1 && AFNI_exist == 0 
    
    D_reorient = vars.ivim.All_D(:,end:-1:1,:,:);
    D_reorient = rot90(D_reorient,1);
    D_resize = imresize(D_reorient,2,'nearest');
    D_zpad = zeros(256,256,26,4);
    D_zpad(:,:,3:24,:) = D_resize;
    D1_zpad = D_zpad(:,:,:,1);
    D2_zpad = D_zpad(:,:,:,2);    
    D3_zpad = D_zpad(:,:,:,3);    
    D4_zpad = D_zpad(:,:,:,4);
    Opt1=struct('Scale',1,'Prefix','D1','Views',[],'verbose',[],'AppendHistory',[],'NoCheck',[],'Slices',[],'Frames',[]);
    WriteBrik(D1_zpad,T1_header,Opt1);
    Opt2=struct('Scale',1,'Prefix','D2a','Views',[],'verbose',[],'AppendHistory',[],'NoCheck',[],'Slices',[],'Frames',[]);
    WriteBrik(D2_zpad,T1_header,Opt2);
    Opt3=struct('Scale',1,'Prefix','D2b','Views',[],'verbose',[],'AppendHistory',[],'NoCheck',[],'Slices',[],'Frames',[]);
    WriteBrik(D3_zpad,T1_header,Opt3);
    Opt4=struct('Scale',1,'Prefix','D3','Views',[],'verbose',[],'AppendHistory',[],'NoCheck',[],'Slices',[],'Frames',[]);
    WriteBrik(D4_zpad,T1_header,Opt4);
    
    
    Dstar_reorient = vars.ivim.All_Dstar(:,end:-1:1,:,:);
    Dstar_reorient = rot90(Dstar_reorient,1);
    Dstar_resize = imresize(Dstar_reorient,2,'nearest');
    Dstar_zpad = zeros(256,256,26,4);
    Dstar_zpad(:,:,3:24,:) = Dstar_resize;
    Dstar1_zpad = Dstar_zpad(:,:,:,1);
    Dstar2_zpad = Dstar_zpad(:,:,:,2);    
    Dstar3_zpad = Dstar_zpad(:,:,:,3);    
    Dstar4_zpad = Dstar_zpad(:,:,:,4);
    Opt1=struct('Scale',1,'Prefix','Dstar1','Views',[],'verbose',[],'AppendHistory',[],'NoCheck',[],'Slices',[],'Frames',[]);
    WriteBrik(Dstar1_zpad,T1_header,Opt1);
    Opt2=struct('Scale',1,'Prefix','Dstar2a','Views',[],'verbose',[],'AppendHistory',[],'NoCheck',[],'Slices',[],'Frames',[]);
    WriteBrik(Dstar2_zpad,T1_header,Opt2);
    Opt3=struct('Scale',1,'Prefix','Dstar2b','Views',[],'verbose',[],'AppendHistory',[],'NoCheck',[],'Slices',[],'Frames',[]);
    WriteBrik(Dstar3_zpad,T1_header,Opt3);
    Opt4=struct('Scale',1,'Prefix','Dstar3','Views',[],'verbose',[],'AppendHistory',[],'NoCheck',[],'Slices',[],'Frames',[]);
    WriteBrik(Dstar4_zpad,T1_header,Opt4);
    
    
    f_reorient = vars.ivim.All_f(:,end:-1:1,:,:);
    f_reorient = rot90(f_reorient,1);
    f_resize = imresize(f_reorient,2,'nearest');
    f_zpad = zeros(256,256,26,4);
    f_zpad(:,:,3:24,:) = f_resize;
    f1_zpad = f_zpad(:,:,:,1);
    f2_zpad = f_zpad(:,:,:,2);    
    f3_zpad = f_zpad(:,:,:,3);    
    f4_zpad = f_zpad(:,:,:,4);
    Opt1=struct('Scale',1,'Prefix','f1','Views',[],'verbose',[],'AppendHistory',[],'NoCheck',[],'Slices',[],'Frames',[]);
    WriteBrik(f1_zpad,T1_header,Opt1);
    Opt2=struct('Scale',1,'Prefix','f2a','Views',[],'verbose',[],'AppendHistory',[],'NoCheck',[],'Slices',[],'Frames',[]);
    WriteBrik(f2_zpad,T1_header,Opt2);
    Opt3=struct('Scale',1,'Prefix','f2b','Views',[],'verbose',[],'AppendHistory',[],'NoCheck',[],'Slices',[],'Frames',[]);
    WriteBrik(f3_zpad,T1_header,Opt3);
    Opt4=struct('Scale',1,'Prefix','f3','Views',[],'verbose',[],'AppendHistory',[],'NoCheck',[],'Slices',[],'Frames',[]);
    WriteBrik(f4_zpad,T1_header,Opt4);
    
    
    
    fDstar_reorient = vars.ivim.All_fDstar(:,end:-1:1,:,:);
    fDstar_reorient = rot90(fDstar_reorient,1);
    fDstar_resize = imresize(fDstar_reorient,2,'nearest');
    fDstar_zpad = zeros(256,256,26,4);
    fDstar_zpad(:,:,3:24,:) = fDstar_resize;
    fDstar1_zpad = fDstar_zpad(:,:,:,1);
    fDstar2_zpad = fDstar_zpad(:,:,:,2);    
    fDstar3_zpad = fDstar_zpad(:,:,:,3);    
    fDstar4_zpad = fDstar_zpad(:,:,:,4);
    Opt1=struct('Scale',1,'Prefix','fDstar1','Views',[],'verbose',[],'AppendHistory',[],'NoCheck',[],'Slices',[],'Frames',[]);
    WriteBrik(fDstar1_zpad,T1_header,Opt1);
    Opt2=struct('Scale',1,'Prefix','fDstar2a','Views',[],'verbose',[],'AppendHistory',[],'NoCheck',[],'Slices',[],'Frames',[]);
    WriteBrik(fDstar2_zpad,T1_header,Opt2);
    Opt3=struct('Scale',1,'Prefix','fDstar2b','Views',[],'verbose',[],'AppendHistory',[],'NoCheck',[],'Slices',[],'Frames',[]);
    WriteBrik(fDstar3_zpad,T1_header,Opt3);
    Opt4=struct('Scale',1,'Prefix','fDstar3','Views',[],'verbose',[],'AppendHistory',[],'NoCheck',[],'Slices',[],'Frames',[]);
    WriteBrik(fDstar4_zpad,T1_header,Opt4);
  
    Opt=struct('Scale',1,'Prefix','T1','Views',[],'verbose',[],'AppendHistory',[],'NoCheck',[],'Slices',[],'Frames',[]);
    WriteBrik(T1,T1_header,Opt);
end