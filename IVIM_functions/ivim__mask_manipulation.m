function ivim__mask_manipulation(iSubj,iPrePost)
disp('Adjusting masks... ')
global vars init_vars

%% Downsample the mask and anatomical data to match resolution of IVIM
mask_og= zeros(128,128,vars.params.slices); mask_og1 = mask_og;
vars.mask.anat = zeros(128,128,vars.params.slices); anat = vars.mask.anat;
MultMask_og= zeros(128,128,vars.params.slices); MultMask_og1 = MultMask_og;
ESMask_og= zeros(128,128,vars.params.slices); ESMask_og1 = ESMask_og;
for j = 1:128
    for k = 1:128
        mask_og(j,k,:) = vars.mask.mask_highres(j*2,k*2,:)/max(vars.mask.mask_highres(:));
        anat(j,k,:) = vars.mask.anat_highres(j*2,k*2,:);
        MultMask_og(j,k,:) = vars.mask.MultMask_highres(j*2,k*2,:)/max(vars.mask.MultMask_highres(:));
        ESMask_og(j,k,:) = vars.mask.ESMask_highres(j*2,k*2,:)/max(vars.mask.ESMask_highres(:));
    end
end
%%
if init_vars.Subj_list{iSubj}=='HC06_01' & iPrePost == 1
    disp('Mask is being adjusted for HC06_01 Pre')
    mask_og1(2:end,1:end-4,:) = mask_og(1:end-1,5:end,:);
    MultMask_og1(2:end,1:end-4,:) = MultMask_og(1:end-1,5:end,:);
    ESMask_og1(2:end,1:end-4,:) = ESMask_og(1:end-1,5:end,:);
elseif init_vars.Subj_list{iSubj}=='HC06_02' & iPrePost == 1
    disp('Mask is being adjusted for HC06_02 Pre')
    mask_og1(1:end-1,8:end,:) = mask_og(2:end,1:end-7,:);
    MultMask_og1(1:end-1,8:end,:) = MultMask_og(2:end,1:end-7,:);
    ESMask_og1(1:end-1,8:end,:) = ESMask_og(2:end,1:end-7,:);
else
    mask_og1(2:end,:,:) = mask_og(1:end-1,:,:);
    MultMask_og1(2:end,:,:) = MultMask_og(1:end-1,:,:);
    ESMask_og1(2:end,:,:) = ESMask_og(1:end-1,:,:);
end

vars.mask.mask_uneroded = mask_og1;
vars.mask.anat = anat;

%% View mask and anatomical data with IVIM data
if vars.saveQAall == 1 || vars.vis == 1
    figure(110)
    imshow3D([squeeze(vars.ivim.b_images_filt3(:,:,:,1)),anat,mask_og1*double(max(vars.ivim.b_images_raw(:))/2),ESMask_og1*double(max(vars.ivim.b_images_raw(:))/2),MultMask_og1*double(max(vars.ivim.b_images_raw(:))/2)],[0 max(vars.ivim.b_images_raw(:))/2]),
    if vars.saveQAall == 1, saveas(gcf,fullfile(vars.dir.QApath,'Presorting','b0_anat_ROIs'),'fig'), end
    if vars.vis == 0, close(figure(110)), end
end

%% Erode Mask
se1 = strel('disk',3);

mask_full = mask_og1;
ESMask_full = ESMask_og1;
MultMask_full = MultMask_og1;

vars.mask.mask_fullrange = imerode(mask_full,se1);
vars.mask.ESMask = ESMask_full.*vars.mask.mask_fullrange;
vars.mask.MultMask = MultMask_full.*vars.mask.mask_fullrange;


%% Variables needed:
maskmin = 1;
maskmax=22;
if init_vars.Subj_list{iSubj}=='HC04_01'
    maskmax=19;
elseif init_vars.Subj_list{iSubj}=='HC04_02'
    maskmax=20;
elseif init_vars.Subj_list{iSubj}=='HC06_02'
    maskmax=19;
elseif init_vars.Subj_list{iSubj}=='HC08_02'
    maskmax=20;
elseif init_vars.Subj_list{iSubj}=='HC12_01'
    maskmax=17;
elseif init_vars.Subj_list{iSubj}=='HC12_02'
    maskmax=17;
elseif init_vars.Subj_list{iSubj}=='HC13_02'
    maskmin = 2;
    maskmax = 19;
elseif init_vars.Subj_list{iSubj}=='HC15_01'
    maskmin = 5;
    maskmax = 22;
    
elseif init_vars.Subj_list{iSubj}=='LS03_01'
    maskmax=20;
elseif init_vars.Subj_list{iSubj}=='LS06_01'
    maskmax=20;
elseif init_vars.Subj_list{iSubj}=='LS14_01'
    maskmax=17;
elseif init_vars.Subj_list{iSubj}=='LS17_01'
    maskmax=21;
elseif init_vars.Subj_list{iSubj}=='LS03_02'
    maskmax=21;
elseif init_vars.Subj_list{iSubj}=='LS14_02'
    maskmin = 7;
    maskmax=19;
end
vars.mask.mask = zeros(vars.params.read,vars.params.pe,vars.params.slices);
vars.mask.mask(:,:,maskmin:maskmax) = vars.mask.mask_fullrange(:,:,maskmin:maskmax);

if init_vars.Subj_list{iSubj}=='HC99_02'
    if iPrePost==1
        vars.mask.mask = vars.mask.mask_fullrange(:,:,maskmax:-1:maskmin);
    end
end


