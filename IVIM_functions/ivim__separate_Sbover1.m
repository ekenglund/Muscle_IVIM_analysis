function ivim__separate_Sbover1(iSubj,iPrePost)

global vars init_vars
%%
for k = 1:vars.params.slices
    for j = 1:length(vars.ivim.b)-1
        SbSo(:,:,k,j) = vars.ivim.Sb(:,:,k,j)./vars.ivim.So(:,:,k);
        SbSo_mask(:,:,k,j) = SbSo(:,:,k,j).*vars.mask.mask(:,:,k);
    end
end


mask_over1 = zeros(128,128,22);
mask_under1 = zeros(128,128,22);

for r = 1:128
    for c = 1:128 
        for s = 1:22
            if SbSo_mask(r,c,s,1)>1
                mask_over1(r,c,s) = 1;
            elseif SbSo_mask(r,c,s,1)>0
                mask_under1(r,c,s) = 1;
            end
        end
    end
end
vars.mask.mask_under1 = mask_under1;
vars.mask.mask_over1 = mask_over1;
% 
% figure, for k = 1:22, imagesc([mask_over1(:,:,k),mask_under1(:,:,k)]), axis image,title(k), pause(0.1),drawnow, end
% %%
size_under1 = length(find(mask_under1));
size_over1 = length(find(mask_over1));
vars.mask.percent_over1 = size_over1/(size_over1+size_under1)*100;

% %%
% for k = 1:4
%     D = vars.ivim.All_D(:,:,:,k);
%     f = vars.ivim.All_f(:,:,:,k);
%     Dstar = vars.ivim.All_Dstar(:,:,:,k);
%     fDstar = vars.ivim.All_fDstar(:,:,:,k);
%     resnorm = vars.ivim.All_resnorm(:,:,:,k);
%  
%     
%     vars.ivim.All_D_under1(:,:,:,k) = D.*mask_under1;
%     vars.ivim.All_f_under1(:,:,:,k) = f.*mask_under1;
%     vars.ivim.All_Dstar_under1(:,:,:,k) = Dstar.*mask_under1;
%     vars.ivim.All_fDstar_under1(:,:,:,k) = fDstar.*mask_under1;
%     vars.ivim.All_resnorm_under1(:,:,:,k) = resnorm.*mask_under1;
%     vars.ivim.All_D_over1(:,:,:,k) = D.*mask_over1;
%     vars.ivim.All_f_over1(:,:,:,k) = f.*mask_over1;
%     vars.ivim.All_Dstar_over1(:,:,:,k) = Dstar.*mask_over1;
%     vars.ivim.All_fDstar_over1(:,:,:,k) = fDstar.*mask_over1;
%     vars.ivim.All_resnorm_over1(:,:,:,k) = resnorm.*mask_over1;
%     
    
    
%     mean_D_under1(k) = mean(D(find(mask_under1)));
%     mean_f_under1(k) = mean(f(find(mask_under1)));
%     mean_Dstar_under1(k) = mean(Dstar(find(mask_under1)));
%     mean_fDstar_under1(k) = mean(fDstar(find(mask_under1)));
%     std_D_under1(k) = std(D(find(mask_under1)));
%     std_f_under1(k) = std(f(find(mask_under1)));
%     std_Dstar_under1(k) = std(Dstar(find(mask_under1)));
%     std_fDstar_under1(k) = std(fDstar(find(mask_under1)));
%     
%     mean_D_over1(k) = mean(D(find(mask_over1)));
%     mean_f_over1(k) = mean(f(find(mask_over1)));
%     mean_Dstar_over1(k) = mean(Dstar(find(mask_over1)));
%     mean_fDstar_over1(k) = mean(fDstar(find(mask_over1)));
%     std_D_over1(k) = std(D(find(mask_over1)));
%     std_f_over1(k) = std(f(find(mask_over1)));
%     std_Dstar_over1(k) = std(Dstar(find(mask_over1)));
%     std_fDstar_over1(k) = std(fDstar(find(mask_over1)));
    
    
%     figure, subplot(2,2,1), hold on, histogram(D(find(mask_under1)),'BinWidth',5e-5)
%     histogram(D(find(mask_over1)),'BinWidth',5e-5), title(['D ',num2str(k)]), 
%     legend(['Under1 = ',num2str(mean_D_under1(k))], ['Over1 = ',num2str(mean_D_over1(k))]),
%     subplot(2,2,2),hold on, histogram(f(find(mask_under1)),'BinWidth',0.01)
%     histogram(f(find(mask_over1)),'BinWidth',0.01), title(['f ',num2str(k)])
%     legend(['Under1 = ',num2str(mean_f_under1(k))], ['Over1 = ',num2str(mean_f_over1(k))]),
%     subplot(2,2,3), hold on, histogram(Dstar(find(mask_under1)),'BinWidth',0.001,'BinLimits',[0 0.025])
%     histogram(Dstar(find(mask_over1)),'BinWidth',0.001,'BinLimits',[0 0.025]), title(['Dstar ',num2str(k)])
%     legend(['Under1 = ',num2str(mean_Dstar_under1(k))], ['Over1 = ',num2str(mean_Dstar_over1(k))]),
%     subplot(2,2,4), hold on, histogram(fDstar(find(mask_under1)),'BinWidth',0.001,'BinLimits',[0 0.02])
%     histogram(fDstar(find(mask_over1)),'BinWidth',0.001,'BinLimits',[0 0.02]), title(['fDstar ',num2str(k)]), 
%     legend(['Under1 = ',num2str(mean_fDstar_under1(k))], ['Over1 = ',num2str(mean_fDstar_over1(k))]),
%     
%     saveas(gcf,fullfile(vars.dir.QApath,'Maps',['Histogram_SeparatedOverUnder_Fit',num2str(k)]),'fig'), 
end
% %%
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% for k = 1:vars.params.slices
%         subplot(2,2,1), imagesc([vars.ivim.All_D_under1(:,:,k,1),vars.ivim.All_D_under1(:,:,k,2),vars.ivim.All_D_under1(:,:,k,3),vars.ivim.All_D_under1(:,:,k,4);...
%             vars.ivim.All_D_over1(:,:,k,1),vars.ivim.All_D_over1(:,:,k,2),vars.ivim.All_D_over1(:,:,k,3),vars.ivim.All_D_over1(:,:,k,4)]), caxis([0 0.0025]), colorbar,axis image, title(['D Slice #',num2str(k)])
%         subplot(2,2,2), imagesc([vars.ivim.All_f_under1(:,:,k,1),vars.ivim.All_f_under1(:,:,k,2),vars.ivim.All_f_under1(:,:,k,3),vars.ivim.All_f_under1(:,:,k,4);...
%             vars.ivim.All_f_over1(:,:,k,1),vars.ivim.All_f_over1(:,:,k,2),vars.ivim.All_f_over1(:,:,k,3),vars.ivim.All_f_over1(:,:,k,4)]), caxis([0 0.2]), colorbar,axis image, title(['f Slice #',num2str(k)])
%         subplot(2,2,3), imagesc([vars.ivim.All_Dstar_under1(:,:,k,1),vars.ivim.All_Dstar_under1(:,:,k,2),vars.ivim.All_Dstar_under1(:,:,k,3),vars.ivim.All_Dstar_under1(:,:,k,4);...
%             vars.ivim.All_Dstar_over1(:,:,k,1),vars.ivim.All_Dstar_over1(:,:,k,2),vars.ivim.All_Dstar_over1(:,:,k,3),vars.ivim.All_Dstar_over1(:,:,k,4)]), caxis([0 0.025]), colorbar,axis image,  title(['Dstar Slice #',num2str(k)])
%         subplot(2,2,4), imagesc([vars.ivim.All_fDstar_under1(:,:,k,1),vars.ivim.All_fDstar_under1(:,:,k,2),vars.ivim.All_fDstar_under1(:,:,k,3),vars.ivim.All_fDstar_under1(:,:,k,4);...
%             vars.ivim.All_fDstar_over1(:,:,k,1),vars.ivim.All_fDstar_over1(:,:,k,2),vars.ivim.All_fDstar_over1(:,:,k,3),vars.ivim.All_fDstar_over1(:,:,k,4)]), caxis([0 0.0025]), colorbar,axis image, title(['fDstar Slice #',num2str(k)])
%        drawnow,
%        saveas(gcf,fullfile(vars.dir.QApath,'Maps','SliceMaps',['SeparatedMaps_TopUnder1_BottomOver1_',num2str(k)]),'fig'),
% end
%     %%
%     
%     
% fit_mean_under(:,1) = mean_D_under1';
% fit_mean_under(:,2) = mean_f_under1';
% fit_mean_under(:,3) = mean_Dstar_under1';
% fit_mean_under(:,4) = mean_fDstar_under1';
% fit_mean_under(1,5) = size_under1;
% 
% fit_mean_over(:,1) = mean_D_over1';
% fit_mean_over(:,2) = mean_f_over1';
% fit_mean_over(:,3) = mean_Dstar_over1';
% fit_mean_over(:,4) = mean_fDstar_over1';
% fit_mean_over(1,5) = size_over1;
% 
% fit_std_under(:,1) = std_D_under1';
% fit_std_under(:,2) = std_f_under1';
% fit_std_under(:,3) = std_Dstar_under1';
% fit_std_under(:,4) = std_fDstar_under1';
% 
% fit_std_over(:,1) = std_D_over1';
% fit_std_over(:,2) = std_f_over1';
% fit_std_over(:,3) = std_Dstar_over1';
% fit_std_over(:,4) = std_fDstar_over1';
% 
% fits_all = zeros(4,23);
% fits_all(:,1:5) = fit_mean_under;
% fits_all(:,7:10) = fit_std_under;
% fits_all(:,13:17) = fit_mean_over;
% fits_all(:,19:22) = fit_std_over;
% fits_all = fits_all';
% 
% % structure is under1 first row 1: D1, D2a, D2b, D3
% % row 2: f1, f2a, f2b, f3
% % row 3: Dstar1, Dstar2a, Dstar2b, Dstar3
% % row 4: fDstar1, fDstar2a, fDstar2b, fDstar3
% % row 5, column 1 = number of voxels under 1
% % then rows 7-10 are standard deviation of under 1
% % then this pattern repeats for mean, standard deviation of voxels with
% % S(b=10)>S0
% 
% xlwrite(fullfile(vars.dir.basedir,'Summary_Results',['20200310_IVIM_ParameterMapsResults_SeparatedOverUnder',num2str(vars.MC_option),'.xls']),fits_all,[init_vars.Subj_list{iSubj},'_',num2str(iPrePost)]);
% 
% % 
% % fit1(:,1) = [D1_mean; f1_mean; Dstar1_mean; fDstar1_mean; r1_mean];
% %     
% %     for k = 1:22, 
% %     figure, imagesc(
% %     xlwrite(fullfile(vars.dir.basedir,'Summary_Results',['20200310_IVIM_ParameterMapsResults_SeparatedOverUnder',num2str(vars.MC_option),'.xls']),fits_mean,[vars.Subj_list{iSubj},'_',num2str(iPrePost)],'A1:E4');
% % 
% % end
