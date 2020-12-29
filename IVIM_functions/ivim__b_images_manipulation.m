function ivim__b_images_manipulation
disp('Sorting b images... ')
global vars


b_images_raw_all_temp = permute(vars.ivim.b_images_raw,[1,2,4,3]);
b_images_raw_all = reshape(b_images_raw_all_temp,[size(vars.ivim.b_images_raw,1), size(vars.ivim.b_images_raw,2), size(vars.ivim.b_images_raw,3)*size(vars.ivim.b_images_raw,4)]);

b_images_filt3_all_temp = permute(vars.ivim.b_images_filt3,[1,2,4,3]);
b_images_filt3_all = reshape(b_images_filt3_all_temp,[size(vars.ivim.b_images_raw,1), size(vars.ivim.b_images_raw,2), size(vars.ivim.b_images_raw,3)*size(vars.ivim.b_images_raw,4)]);

if vars.saveQAall == 1 || vars.vis == 1
    figure(120), imshow3D([b_images_raw_all,b_images_filt3_all],[0 max(vars.ivim.b_images_raw(:))/5])
    if vars.saveQAall == 1, saveas(gcf,fullfile(vars.dir.QApath,'Presorting','b_images_raw_filt3'),'fig'), end
    if vars.vis == 0, close(figure(120)), end
end
clear b_images_raw_all_temp b_images_raw_all b_images_filt3_all_temp b_images_filt3_all b_images_reg_all_temp b_images_reg_all

%% sort images based on tensor file
[vars.ivim.Sb, vars.ivim.So, vars.ivim.b_images, vars.ivim.b_0, vars.ivim.b_images_all,vars.ivim.b_images_time, vars.ivim.b, vars.ivim.b_time] = tensor_sorter(vars.ivim.b_images_filt3, vars.params.tensor, vars.params.dirs);
%for 149 directions, b_0(:,:,:,:,2:5) correspond temporally to
%the b_images data
%Sb is the averaged data over repeats
%So is the averaged non-diffusion-weighted data
%b_images is the unaveraged diffusion-weighted data
%b_0 is the unaveraged non-diffusion-weighted data
%b_images_time is all data (averaged xyz diffusion) over time
%b is the vector of b-values

if vars.saveQAall == 1 || vars.vis == 1
    figure(121), for k = 1:vars.params.slices
        for j = 1:size(vars.ivim.b_images,4)
            imagesc(vars.ivim.Sb(:,:,k,j)), hold on, axis image, title([k,j]), contour(vars.mask.mask(:,:,k),'k'), caxis([0 max(vars.ivim.b_images(:))/4]), hold off,  drawnow;
        end
    end
    
    Sb_all_temp = permute(vars.ivim.Sb,[1,2,4,3]);
    Sb_all = reshape(Sb_all_temp,[size(vars.ivim.Sb,1), size(vars.ivim.Sb,2), size(vars.ivim.Sb,3)*size(vars.ivim.Sb,4)]);
    figure(122), imshow3D(Sb_all,[0 max(vars.ivim.b_images_raw(:))/5])
    if vars.saveQAall == 1, saveas(gcf,fullfile(vars.dir.QApath,'Maps','Sb_images'),'fig'), end
    if vars.vis == 0, close(figure(122)), end
    clear Sb_all_temp Sb_all
end