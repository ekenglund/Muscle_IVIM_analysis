function ivim__fit_resolve_time_varymask(iSubj,iPrePost)
disp('Computing IVIM parameters in full ROI as a function of time... ')
global vars init_vars

for imask = 1:4
    % apply mask to full muscle (e.g. not resolved slice-wise)
    for k = 1:size(vars.ivim.b_images_time,4)
        b_images_temp = double(squeeze(vars.ivim.b_images_time(:,:,:,k)));
        if imask == 1
            mask_temp = vars.mask.fullexclusion.mask1;
        elseif imask == 2
            mask_temp = vars.mask.fullexclusion.mask2a;
        elseif imask == 3
            mask_temp = vars.mask.fullexclusion.mask2b;
        elseif imask == 4
            mask_temp = vars.mask.fullexclusion.mask3;
        end
        Sb_bytime(k) = mean(b_images_temp(find(mask_temp)),'omitnan');
        %     Sb_ES(k) = mean(b_images_temp(find(vars.mask.ESMask)),'omitnan');
        %     Sb_Mult(k) = mean(b_images_temp(find(vars.mask.MultMask)),'omitnan');
    end
    
    % if vars.saveQA == 1 || vars.vis == 1
    %     figure(140), hold on,grid on, for k = 1:size(vars.ivim.b_images_all,5)
    %         plot(vars.ivim.b,Sb_bytime(((k-1)*13+1:k*13)),'.-','Color',[0 0 2/(k+1)])
    %     end
    %     for k = 1:size(vars.ivim.b_images_all,5)
    %         plot(vars.ivim.b,Sb_ES(((k-1)*13+1:k*13)),'.-','Color',[0 2/(k+1) 0])
    %     end
    %     for k = 1:size(vars.ivim.b_images_all,5)
    %         plot(vars.ivim.b,Sb_Mult(((k-1)*13+1:k*13)),'.-','Color',[2/(k+1) 0 0])
    %     end
    %     legend('PS1','PS2','PS3','PS4','ES1','ES2','ES3','ES4','Mult1','Mult2','Mult3','Mult4')
    %     title('S(b) averaged over entire muscle ROI')
    %     xlabel('b (s/mm^2)'), ylabel('Signal intensity (AU)')
    %     if vars.saveQA == 1, saveas(gcf,fullfile(vars.dir.QApath,'FitMean','FullROI_timeresolved_Sb'),'fig'), end
    %     if vars.vis == 0, close(figure(140)), end
    % end
    
    SbSo_bytime = zeros(size(vars.ivim.b_images_all,5),size(vars.ivim.b,2)-1);
    % SbSo_ES = SbSo_bytime;
    % SbSo_Mult = SbSo_bytime;
    for k = 1:size(vars.ivim.b_images_all,5)
        SbSo_bytime(k,:) = Sb_bytime((k-1)*13+2:k*13)./(Sb_bytime((k-1)*13+1));
        %     SbSo_ES(k,:) = Sb_ES((k-1)*13+2:k*13)./(Sb_ES((k-1)*13+1));
        %     SbSo_Mult(k,:) = Sb_Mult((k-1)*13+2:k*13)./(Sb_Mult((k-1)*13+1));
    end
    
    % if vars.saveQA == 1 || vars.vis == 1
    %     figure(141), hold on,grid on, for k = 1:size(vars.ivim.b_images_all,5)
    %         plot(vars.ivim.b,[1 squeeze(SbSo_bytime(k,:))],'.-','Color',[0 0 2/(k+1)])
    %     end
    %     for k = 1:size(vars.ivim.b_images_all,5)
    %         plot(vars.ivim.b,[1 squeeze(SbSo_ES(k,:))],'.-','Color',[0 2/(k+1) 0])
    %     end
    %     for k = 1:size(vars.ivim.b_images_all,5)
    %         plot(vars.ivim.b,[1 squeeze(SbSo_Mult(k,:))],'.-','Color',[2/(k+1) 0 0])
    %     end
    %     legend('PS1','PS2','PS3','PS4','ES1','ES2','ES3','ES4','Mult1','Mult2','Mult3','Mult4')
    %     title('S(b)/So averaged over entire muscle ROI')
    %     xlabel('b (s/mm^2)'), ylabel('Normalized signal intensity (AU)')
    %     if vars.saveQA == 1, saveas(gcf,fullfile(vars.dir.QApath,'FitMean','Timewise_fullROIaverage_SbSo'),'fig'), end
    %     if vars.vis == 0, close(figure(141)), end
    % end
    
    for k = 1:size(vars.ivim.b_images_all,5)
        [coeffs(:,:,k), resnorm(:,k)] = ivim__parameter_fitting(SbSo_bytime(k,:));
    end
    
    vars.ivim.bytime.SbSo_bytime = SbSo_bytime;
    % matrix is [f-t1    f-t2    f-t3    f-t4
    %            D*-t1   D*-t2   D*-t3   D*-t4
    %            D-t1    D-t2    D-t3    D-t4
    %            fD*-t1  fD*-t2  fD*-t3  fD*-t4
    %            r-t1    r-t2    r-t3    r-t4 ] ** r is resnorm
    if imask == 1
        vars.ivim.bytime.coeffs_1partfit_fullexclusion = squeeze(coeffs(1,:,:));
        vars.ivim.bytime.coeffs_1partfit_fullexclusion(4,:) = vars.ivim.bytime.coeffs_1partfit_fullexclusion(1,:).*vars.ivim.bytime.coeffs_1partfit_fullexclusion(2,:);
        vars.ivim.bytime.coeffs_1partfit_fullexclusion(5,:) = resnorm(1,:);
    elseif imask == 2
        vars.ivim.bytime.coeffs_2partfita_fullexclusion = squeeze(coeffs(2,:,:));
        vars.ivim.bytime.coeffs_2partfita_fullexclusion(4,:) = vars.ivim.bytime.coeffs_2partfita_fullexclusion(1,:).*vars.ivim.bytime.coeffs_2partfita_fullexclusion(2,:);
        vars.ivim.bytime.coeffs_2partfita_fullexclusion(5,:) = resnorm(2,:);
    elseif imask == 3
        vars.ivim.bytime.coeffs_2partfitb_fullexclusion = squeeze(coeffs(3,:,:));
        vars.ivim.bytime.coeffs_2partfitb_fullexclusion(4,:) = vars.ivim.bytime.coeffs_2partfitb_fullexclusion(1,:).*vars.ivim.bytime.coeffs_2partfitb_fullexclusion(2,:);
        vars.ivim.bytime.coeffs_2partfitb_fullexclusion(5,:) = resnorm(3,:);
    elseif imask == 4
        vars.ivim.bytime.coeffs_3partfit_fullexclusion = squeeze(coeffs(4,:,:));
        vars.ivim.bytime.coeffs_3partfit_fullexclusion(4,:) = vars.ivim.bytime.coeffs_3partfit_fullexclusion(1,:).*vars.ivim.bytime.coeffs_3partfit_fullexclusion(2,:);
        vars.ivim.bytime.coeffs_3partfit_fullexclusion(5,:) = resnorm(4,:);
    end
end

coeffs_1partfit = vars.ivim.bytime.coeffs_1partfit_fullexclusion';
coeffs_2partfita = vars.ivim.bytime.coeffs_2partfita_fullexclusion';
coeffs_2partfitb = vars.ivim.bytime.coeffs_2partfitb_fullexclusion';
coeffs_3partfit = vars.ivim.bytime.coeffs_3partfit_fullexclusion';

% %%
% if vars.saveQA == 1 || vars.vis == 1
%     coeffs_1partfit = vars.ivim.bytime.coeffs_1partfit';
%     coeffs_2partfita = vars.ivim.bytime.coeffs_2partfita';
%     coeffs_2partfitb = vars.ivim.bytime.coeffs_2partfitb';
%     coeffs_3partfit = vars.ivim.bytime.coeffs_3partfit';
%     b = vars.ivim.b;
%     for k = 1:4
%         SbSo_fit1 = coeffs_1partfit(k,1)*exp(-b*coeffs_1partfit(k,2))+(1-coeffs_1partfit(k,1))*exp(-b*coeffs_1partfit(k,3));
%         SbSo_fit2 = coeffs_2partfita(k,1)*exp(-b*coeffs_2partfita(k,2))+(1-coeffs_2partfita(k,1))*exp(-b*coeffs_2partfita(k,3));
%         SbSo_fit3 = coeffs_2partfitb(k,1)*exp(-b*coeffs_2partfitb(k,2))+(1-coeffs_2partfitb(k,1))*exp(-b*coeffs_2partfitb(k,3));
%         SbSo_fit4 = coeffs_3partfit(k,1)*exp(-b*coeffs_3partfit(k,2))+(1-coeffs_3partfit(k,1))*exp(-b*coeffs_3partfit(k,3));
%         figure(142), subplot(1,4,k), plot(b,[1 SbSo_bytime(k,:)],'ko-','LineWidth',1.5,'MarkerSize',8), grid on, hold on
%         plot(b,SbSo_fit1,'.-','MarkerSize',10), plot(b,SbSo_fit2,'.-','MarkerSize',10), plot(b,SbSo_fit3,'.-','MarkerSize',10), plot(b,SbSo_fit4,'.-','MarkerSize',10)
%         title(['Repeat',num2str(k)]), axis([0 700 0 1]), if k == 1, ylabel('Paraspinal ROI'), end
%         legend('Measured','fit1','fit2a','fit2b','fit3')
%         if vars.saveQA == 1, saveas(gcf,fullfile(vars.dir.QApath,'FitMean','Timewise_fullROIaverage_Fits'),'fig'),end
%     end
%     if vars.vis == 0, close(figure(142)), end
% end

%%
if vars.saveQA == 1 || vars.vis == 1
    figure(143), for k = 1:5, subplot(1,5,k),hold on, plot(squeeze(coeffs_1partfit(:,k)),'.-'), grid on,
        if k == 1, title('f'), elseif k == 2, title('D*'), elseif k == 3, title('D'),...
        elseif k == 4, title('fD*'), elseif k == 5, title('resnorm'), end, end
for k = 1:5, subplot(1,5,k),hold on,plot(squeeze(coeffs_2partfita(:,k)),'.-'), grid on, end
for k = 1:5, subplot(1,5,k),hold on, plot(squeeze(coeffs_2partfitb(:,k)),'.-'), grid on, end
for k = 1:5, subplot(1,5,k),hold on, plot(squeeze(coeffs_3partfit(:,k)),'.-'), grid on,
    if k == 1, ylabel('Paraspinal ROI'), end, end
legend('1','2a','2b','3')
if vars.saveQA == 1, saveas(gcf,fullfile(vars.dir.QApath,'FitMean','Timewise_fullROIaverage_IVIMCoeffs_vs_time_ExMask_rmAB'),'fig'), end
if vars.vis == 0, close(figure(143)), end
end


%%
fits = zeros(23,4);
% rows = f, D*, D, fD*, resnorm
% cols = timepoint (1, 2, 3, 4)
% sections = fits 1, 2a, 2b, 3

fits(1:5,:) = vars.ivim.bytime.coeffs_1partfit_fullexclusion;
fits(7:11,:) = vars.ivim.bytime.coeffs_2partfita_fullexclusion;
fits(13:17,:) = vars.ivim.bytime.coeffs_2partfitb_fullexclusion;
fits(19:23,:) = vars.ivim.bytime.coeffs_3partfit_fullexclusion;

% xlwrite(fullfile(vars.dir.basedir,'Summary_Results',['20200721_IVIM_ResolveTimeResults_ExMask_rmAB_',num2str(vars.MC_option),'.xls']),fits,[init_vars.Subj_list{iSubj},'_',num2str(iPrePost)]);


