function ivim__relativeT2
disp('Computing relative T2 changes over time... ')
global vars

% mask b_0
for k = 1:vars.params.slices
    for j = 1:size(vars.ivim.b_0,4)
        mask_temp=vars.mask.mask(:,:,k);
        b_0_temp = double(squeeze(vars.ivim.b_0(:,:,k,j)));
        b_0_mask(:,:,k,j)=b_0_temp.*mask_temp;
        mask_temp = zeros(vars.params.read,vars.params.pe);
    end
end
clear mask_temp b_0_temp
%

for k = 1:vars.params.slices, for j = 1:size(b_0_mask,4)-1
        rel_b0(:,:,k,j) = imdivide(b_0_mask(:,:,k,j+1),b_0_mask(:,:,k,1));
    end
end
rel_b0_all_temp = permute(rel_b0,[1,2,4,3]);
vars.t2.rel_b0_all = reshape(rel_b0_all_temp,[size(rel_b0,1), size(rel_b0,2), size(rel_b0,3)*size(rel_b0,4)]);

figure(130), imshow3D(vars.t2.rel_b0_all,[0.8,1.2])
if vars.saveQAall == 1, saveas(gcf,fullfile(vars.dir.QApath,'T2','rel_b0'),'fig'), end
if vars.vis == 0, close(figure(130)), end
%
% apply mask
for k = 1:vars.params.slices
    for j = 1:size(rel_b0,4)
        mask_temp=vars.mask.mask(:,:,k);
        ESmask_temp = vars.mask.ESMask(:,:,k);
        Mmask_temp = vars.mask.MultMask(:,:,k);
        rel_b0_temp = double(squeeze(rel_b0(:,:,k,j)));
        vars.t2.rel_b0_img_PS(:,:,k,j)=rel_b0_temp.*mask_temp;
        vars.t2.rel_b0_PS_byslice(k,j) = mean(rel_b0_temp(find(mask_temp)));
        vars.t2.rel_b0_img_ES(:,:,k,j)=rel_b0_temp.*ESmask_temp;
        vars.t2.rel_b0_ES_byslice(k,j) = mean(rel_b0_temp(find(ESmask_temp)));
        vars.t2.rel_b0_img_Mult(:,:,k,j)=rel_b0_temp.*Mmask_temp;
        vars.t2.rel_b0_Mult_byslice(k,j) = mean(rel_b0_temp(find(Mmask_temp)));
        mask_temp = zeros(vars.params.read,vars.params.pe); ESmask_temp = mask_temp; Mmask_temp = mask_temp;
    end
end
clear mask_temp rel_b0_temp

% apply mask to full muscle (e.g. not resolved slice-wise)
for k = 1:size(rel_b0,4)
    rel_b0_images_temp = double(squeeze(rel_b0(:,:,:,k)));
    vars.t2.rel_b0_PS_med(k) = median(rel_b0_images_temp(find(vars.mask.mask)),'omitnan');
    vars.t2.rel_b0_ES_med(k) = median(rel_b0_images_temp(find(vars.mask.ESMask)),'omitnan');
    vars.t2.rel_b0_Mult_med(k) = median(rel_b0_images_temp(find(vars.mask.MultMask)),'omitnan');
    vars.t2.rel_b0_PS(k) = mean(rel_b0_images_temp(find(vars.mask.mask)),'omitnan');
    vars.t2.rel_b0_ES(k) = mean(rel_b0_images_temp(find(vars.mask.ESMask)),'omitnan');
    vars.t2.rel_b0_Mult(k) = mean(rel_b0_images_temp(find(vars.mask.MultMask)),'omitnan');
end
for k = 1:size(b_0_mask,4)
    b_0_mask_temp = double(squeeze(b_0_mask(:,:,:,k)));
    vars.t2.b0_PS(k) = mean(b_0_mask_temp(find(vars.mask.mask)),'omitnan');
    vars.t2.b0_ES(k) = mean(b_0_mask_temp(find(vars.mask.ESMask)),'omitnan');
    vars.t2.b0_Mult(k) = mean(b_0_mask_temp(find(vars.mask.MultMask)),'omitnan');
end
if vars.saveQAall == 1 || vars.vis == 1
    figure(131), hold on, grid on
    plot(vars.t2.rel_b0_PS,'b.-')
    plot(vars.t2.rel_b0_ES,'g.-')
    plot(vars.t2.rel_b0_Mult,'r.-')
    title('Relative T2 changes over time')
    ylabel('Normalized signal intensity'), xlabel('Time point')
    legend('PS','ES','Mult')
    if vars.saveQAall == 1, saveas(gcf,fullfile(vars.dir.QApath,'T2','rel_b0_mean'),'fig'), end
    if vars.vis == 0, close(figure(131)), end
end

if vars.saveQAall == 1 || vars.vis == 1
    figure(132), hold on,
    subplot(1,3,1), hold on, for k = 1:vars.params.slices, plot(vars.t2.rel_b0_PS_byslice(k,:),'.-','Color',[0/255 204/(204+k*30) 55/(55+k)]), end, title('DeltaT2 by slice, PS'), ylabel('Normalized SI'), xlabel('Time point'), grid on
    subplot(1,3,2), hold on, for k = 1:vars.params.slices, plot(vars.t2.rel_b0_ES_byslice(k,:),'.-','Color',[0/255 204/(204+(k-1)*20) 0/255]), end, title('DeltaT2 by slice, ES'), ylabel('Normalized SI'), xlabel('Time point'), grid on
    subplot(1,3,3), hold on, for k = 1:vars.params.slices, plot(vars.t2.rel_b0_Mult_byslice(k,:),'.-','Color',[255/(255+(k-1)*30) 50/(255+(k-1)*100) 0/255]), end, title('DeltaT2 by slice, Mult'), ylabel('Normalized SI'), xlabel('Time point'), grid on
    
    % here lighter colors signify lower slice number (more
    % superior)
    
    if vars.saveQAall == 1, saveas(gcf,fullfile(vars.dir.QApath,'T2','rel_b0_mean_byslice'),'fig'), end
    if vars.vis == 0, close(figure(132)), end
end