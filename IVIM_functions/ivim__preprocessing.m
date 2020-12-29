function ivim__preprocessing
disp('Preprocessing... ')

global vars

if vars.saveQAall == 1 || vars.vis == 1
    figure(100), imshow3D(squeeze(vars.ivim.b_images_raw(:,:,:,1)),[0 max(vars.ivim.b_images_raw(:))/2])
    if vars.saveQAall == 1, saveas(gcf,fullfile(vars.dir.QApath,'Presorting','b0_raw'),'fig'), end, 
    if vars.vis == 0, close(figure(100)), end
end

for j = 1:vars.params.slices
    for k = 1:vars.params.dirs+1
        vars.ivim.b_images_filt3(:,:,j,k) = medfilt2(vars.ivim.b_images_raw(:,:,j,k),[3 3]);
    end
end
if vars.saveQAall == 1 || vars.vis == 1
    figure(101), imshow3D(squeeze(vars.ivim.b_images_filt3(:,:,:,6)),[0 max(vars.ivim.b_images_raw(:))/2]),
    if vars.saveQAall == 1, saveas(gcf,fullfile(vars.dir.QApath,'Presorting','b0_filt3'),'fig'), end
    if vars.vis == 0, close(figure(101)), end
end