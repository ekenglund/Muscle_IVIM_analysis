function ivim__fit_resolve_slice_varymask(iSubj,iPrePost)
disp('Computing IVIM parameters in each slice... ')
global vars init_vars
%%
b = vars.ivim.b;
for imask = 1:4
    % apply mask slice-wise (e.g. not resolved by time)
    for k = 1:size(vars.ivim.Sb,3)
        %     mask_temp = vars.mask.mask(:,:,k);
        if imask == 1
            mask_temp = vars.mask.fullexclusion.mask1(:,:,k);
        elseif imask == 2
            mask_temp = vars.mask.fullexclusion.mask2a(:,:,k);
        elseif imask == 3
            mask_temp = vars.mask.fullexclusion.mask2b(:,:,k);
        elseif imask == 4
            mask_temp = vars.mask.fullexclusion.mask3(:,:,k);
        end
            
        So_temp = vars.ivim.So(:,:,k);
        So(k,1) = mean(So_temp(find(mask_temp)),'omitnan');
        for j = 1:size(vars.ivim.Sb,4)
            Sb_temp = vars.ivim.Sb(:,:,k,j);
            Sb(k,j) = mean(Sb_temp(find(mask_temp)),'omitnan');
        end
    end
    S_all_byslice(:,1) = So;
    S_all_byslice(:,2:13) = Sb;
    
    % if vars.saveQA == 1 || vars.vis == 1
    %     figure(150), hold on, grid on, plot(b,S_all_byslice')
    %     xlabel('b (s/mm^2)'), ylabel('Signal Intensity (AU)')
    %     title('Slice-wise mean of S(b) over b')
    %     if vars.saveQA == 1, saveas(gcf,fullfile(vars.dir.QApath,'FitMean','SlicewiseROI_Sb_vs_b'),'fig'),end
    %     if vars.vis == 0, close(figure(150)), end
    % end
    
    SbSo_byslice = zeros(vars.params.slices,length(b)-1);
    for k = 1:length(b)-1
        SbSo_byslice(:,k) = Sb(:,k)./So;
    end
    
    % if vars.saveQA == 1 || vars.vis == 1
    %     figure(151), hold on,grid on, for k = 1:vars.params.slices
    %         plot(vars.ivim.b,[1 squeeze(SbSo_byslice(k,:))],'.-','Color',[0 0 2/(k+1)])
    %     end
    %     title('S(b)/So averaged over time by slice')
    %     xlabel('b (s/mm^2)'), ylabel('Normalized signal intensity (AU)')
    %     if vars.saveQA == 1, saveas(gcf,fullfile(vars.dir.QApath,'FitMean','SlicewiseROI_timeaveraged_SbSo'),'fig'), end
    %     if vars.vis == 0, close(figure(151)), end
    % end
    
    %%
    for k = 1:vars.params.slices
        [coeffs(:,:,k), resnorm(:,k)] = ivim__parameter_fitting(SbSo_byslice(k,:));
    end
    
    vars.ivim.byslice.SbSo_byslice_All = ones(22,13);
    vars.ivim.byslice.SbSo_byslice_All(:,2:13) = SbSo_byslice;
    
    vars.ivim.byslice.SbSo_byslice = SbSo_byslice;
    
    % matrix is [f-t1    f-t2    f-t3    f-t4
    %            D*-t1   D*-t2   D*-t3   D*-t4
    %            D-t1    D-t2    D-t3    D-t4
    %            fD*-t1  fD*-t2  fD*-t3  fD*-t4
    %            r-t1    r-t2    r-t3    r-t4 ] ** r is resnorm
    if imask == 1
        vars.ivim.byslice.coeffs_1partfit_fullexclusion = squeeze(coeffs(1,:,:));
        vars.ivim.byslice.coeffs_1partfit_fullexclusion(4,:) = vars.ivim.byslice.coeffs_1partfit_fullexclusion(1,:).*vars.ivim.byslice.coeffs_1partfit_fullexclusion(2,:);
        vars.ivim.byslice.coeffs_1partfit_fullexclusion(5,:) = resnorm(1,:);
    elseif imask == 2
        vars.ivim.byslice.coeffs_2partfita_fullexclusion = squeeze(coeffs(2,:,:));
        vars.ivim.byslice.coeffs_2partfita_fullexclusion(4,:) = vars.ivim.byslice.coeffs_2partfita_fullexclusion(1,:).*vars.ivim.byslice.coeffs_2partfita_fullexclusion(2,:);
        vars.ivim.byslice.coeffs_2partfita_fullexclusion(5,:) = resnorm(2,:);
    elseif imask == 3
        vars.ivim.byslice.coeffs_2partfitb_fullexclusion = squeeze(coeffs(3,:,:));
        vars.ivim.byslice.coeffs_2partfitb_fullexclusion(4,:) = vars.ivim.byslice.coeffs_2partfitb_fullexclusion(1,:).*vars.ivim.byslice.coeffs_2partfitb_fullexclusion(2,:);
        vars.ivim.byslice.coeffs_2partfitb_fullexclusion(5,:) = resnorm(3,:);
    elseif imask == 4
        vars.ivim.byslice.coeffs_3partfit_fullexclusion = squeeze(coeffs(4,:,:));
        vars.ivim.byslice.coeffs_3partfit_fullexclusion(4,:) = vars.ivim.byslice.coeffs_3partfit_fullexclusion(1,:).*vars.ivim.byslice.coeffs_3partfit_fullexclusion(2,:);
        vars.ivim.byslice.coeffs_3partfit_fullexclusion(5,:) = resnorm(4,:);
    end
end
% 
% 
% 
    coeffs_1partfit = vars.ivim.byslice.coeffs_1partfit_fullexclusion';
    coeffs_2partfita = vars.ivim.byslice.coeffs_2partfita_fullexclusion';
    coeffs_2partfitb = vars.ivim.byslice.coeffs_2partfitb_fullexclusion';
    coeffs_3partfit = vars.ivim.byslice.coeffs_3partfit_fullexclusion';
%%
% if vars.saveQA == 1 || vars.vis == 1

%     b = vars.ivim.b;
%     for k = 1:22
%         SbSo_fit1 = coeffs_1partfit(k,1)*exp(-b*coeffs_1partfit(k,2))+(1-coeffs_1partfit(k,1))*exp(-b*coeffs_1partfit(k,3));
%         SbSo_fit2 = coeffs_2partfita(k,1)*exp(-b*coeffs_2partfita(k,2))+(1-coeffs_2partfita(k,1))*exp(-b*coeffs_2partfita(k,3));
%         SbSo_fit3 = coeffs_2partfitb(k,1)*exp(-b*coeffs_2partfitb(k,2))+(1-coeffs_2partfitb(k,1))*exp(-b*coeffs_2partfitb(k,3));
%         SbSo_fit4 = coeffs_3partfit(k,1)*exp(-b*coeffs_3partfit(k,2))+(1-coeffs_3partfit(k,1))*exp(-b*coeffs_3partfit(k,3));
%         figure(152), subplot(4,6,k), plot(b,[1 SbSo_byslice(k,:)],'ko-','LineWidth',1.5,'MarkerSize',8), grid on, hold on
%         plot(b,SbSo_fit1,'.-','MarkerSize',10), plot(b,SbSo_fit2,'.-','MarkerSize',10), plot(b,SbSo_fit3,'.-','MarkerSize',10), plot(b,SbSo_fit4,'.-','MarkerSize',10)
%         axis([0 700 0 1])
%         title(['Slice',num2str(k)])
%         if k == 22
%             legend('Measured','fit1','fit2a','fit2b','fit3')
%         end
%         if vars.saveQA == 1, saveas(gcf,fullfile(vars.dir.QApath,'FitMean','SlicewiseROI_Fits_ExMask_rmAB'),'fig'),end
%     end
%     if vars.vis == 0, close(figure(152)), end
% end

%%
if vars.saveQA == 1 || vars.vis == 1
    figure(153), for k = 1:5, subplot(1,5,k),hold on, plot(squeeze(coeffs_1partfit(:,k)),'.-'), grid on,
        if k == 1, title('f'), elseif k == 2, title('D*'), elseif k == 3, title('D'),...
        elseif k == 4, title('fD*'), elseif k == 5, title('resnorm'), end, end
for k = 1:5, subplot(1,5,k),hold on,plot(squeeze(coeffs_2partfita(:,k)),'.-'), grid on, end
for k = 1:5, subplot(1,5,k),hold on, plot(squeeze(coeffs_2partfitb(:,k)),'.-'), grid on, end
for k = 1:5, subplot(1,5,k),hold on, plot(squeeze(coeffs_3partfit(:,k)),'.-'), grid on,
    if k == 1, ylabel('Paraspinal ROI'), end, end
legend('1','2a','2b','3')
if vars.saveQA == 1, saveas(gcf,fullfile(vars.dir.QApath,'FitMean','SlicewiseROI_IVIMCoeffs_vs_slice_ExMask_rmAB'),'fig'), end
if vars.vis == 0, close(figure(153)), end
end
% 
% 
% %%
% fits = zeros(23,22); 
% % rows = f, D*, D, fD*, resnorm 
% % cols = timepoint (1, 2, 3, 4)
% % sections = fits 1, 2a, 2b, 3
% 
fits(1:5,:) = vars.ivim.byslice.coeffs_1partfit_fullexclusion;
fits(7:11,:) = vars.ivim.byslice.coeffs_2partfita_fullexclusion;
fits(13:17,:) = vars.ivim.byslice.coeffs_2partfitb_fullexclusion;
fits(19:23,:) = vars.ivim.byslice.coeffs_3partfit_fullexclusion;
% 
% xlwrite(fullfile(vars.dir.basedir,'Summary_Results',['20200721_IVIM_ResolveSliceResults_ExMask_rmAB_',num2str(vars.MC_option),'.xls']),fits,[init_vars.Subj_list{iSubj},'_',num2str(iPrePost)]);
% 
% 
% 
% 
