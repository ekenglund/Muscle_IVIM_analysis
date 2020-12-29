function ivim__fit_map(iSubj,iPrePost)
disp('Computing IVIM parameter maps... ')
%%
global vars init_vars

b = vars.ivim.b;

Sb = vars.ivim.Sb;
So = vars.ivim.So;

mask = vars.mask.mask_under1;
% se1 = strel('disk',3);
% dilmask = imdilate(vars.mask.mask_og1,se1);

for k = 1:vars.params.slices
    for j = 1:length(b)-1
        SbSo(:,:,k,j) = imdivide(Sb(:,:,k,j),So(:,:,k));
        SbSo_masked(:,:,k,j) = SbSo(:,:,k,j).*mask(:,:,k);
    end
end

D1 = zeros(vars.params.read,vars.params.pe,vars.params.slices); D2a = D1; D2b = D1; D3 = D1;
f1 = zeros(vars.params.read,vars.params.pe,vars.params.slices); f2a = D1; f2b = D1; f3 = D1;
Dstar1 = zeros(vars.params.read,vars.params.pe,vars.params.slices); Dstar2a = D1; Dstar2b = D1; Dstar3 = D1;
fDstar1 = zeros(vars.params.read,vars.params.pe,vars.params.slices); fDstar2a = D1; fDstar2b = D1; fDstar3 = D1;
r1 = zeros(vars.params.read,vars.params.pe,vars.params.slices); r2a = D1; r2b = D1; r3 = D1; %resnorm

%%
idx = find(mask);
h = waitbar(0,['Analyzing IVIM Maps for ',init_vars.Subj_list(iSubj),num2str(iPrePost)]);
for k = 1:length(idx)
    [r,c,s] = ind2sub(size(mask),idx(k));
    [coeffs, resnorm] = ivim__parameter_fitting(squeeze(SbSo(r,c,s,:))');
    f1(r,c,s) = coeffs(1,1); Dstar1(r,c,s) = coeffs(1,2); D1(r,c,s) = coeffs(1,3); fDstar1(r,c,s) = coeffs(1,1).*coeffs(1,2); r1(r,c,s) = resnorm(1);
    f2a(r,c,s) = coeffs(2,1); Dstar2a(r,c,s) = coeffs(2,2); D2a(r,c,s) = coeffs(2,3); fDstar2a(r,c,s) = coeffs(2,1).*coeffs(2,2); r2a(r,c,s) = resnorm(2);
    f2b(r,c,s) = coeffs(3,1); Dstar2b(r,c,s) = coeffs(3,2); D2b(r,c,s) = coeffs(3,3); fDstar2b(r,c,s) = coeffs(3,1).*coeffs(3,2); r2b(r,c,s) = resnorm(3);
    f3(r,c,s) = coeffs(4,1); Dstar3(r,c,s) = coeffs(4,2); D3(r,c,s) = coeffs(4,3); fDstar3(r,c,s) = coeffs(4,1).*coeffs(4,2); r3(r,c,s) = resnorm(4);
    waitbar(k/length(idx),h)
end
%%

vars.ivim.All_D = zeros(vars.params.read,vars.params.pe,vars.params.slices,4);
vars.ivim.All_f = zeros(vars.params.read,vars.params.pe,vars.params.slices,4);
vars.ivim.All_Dstar = zeros(vars.params.read,vars.params.pe,vars.params.slices,4);
vars.ivim.All_fDstar = zeros(vars.params.read,vars.params.pe,vars.params.slices,4);
vars.ivim.All_resnorm = zeros(vars.params.read,vars.params.pe,vars.params.slices,4);

vars.ivim.All_D(:,:,:,1) = D1; vars.ivim.All_D(:,:,:,2) = D2a; vars.ivim.All_D(:,:,:,3) = D2b; vars.ivim.All_D(:,:,:,4) = D3;
vars.ivim.All_f(:,:,:,1) = f1; vars.ivim.All_f(:,:,:,2) = f2a; vars.ivim.All_f(:,:,:,3) = f2b; vars.ivim.All_f(:,:,:,4) = f3;
vars.ivim.All_Dstar(:,:,:,1) = Dstar1; vars.ivim.All_Dstar(:,:,:,2) = Dstar2a; vars.ivim.All_Dstar(:,:,:,3) = Dstar2b; vars.ivim.All_Dstar(:,:,:,4) = Dstar3;
vars.ivim.All_fDstar(:,:,:,1) = fDstar1; vars.ivim.All_fDstar(:,:,:,2) = fDstar2a; vars.ivim.All_fDstar(:,:,:,3) = fDstar2b; vars.ivim.All_fDstar(:,:,:,4) = fDstar3;
vars.ivim.All_resnorm(:,:,:,1) = r1; vars.ivim.All_resnorm(:,:,:,2) = r2a; vars.ivim.All_resnorm(:,:,:,3) = r2b; vars.ivim.All_resnorm(:,:,:,4) = r3;

if vars.saveQA == 1 || vars.vis == 1
    figure(160)
    imshow3D([D1,D2a,D2b,D3],[0 2.5e-3]), colorbar, title('D'), if vars.saveQA==1, saveas(gcf,fullfile(vars.dir.QApath,'Maps','D'),'fig'), end
    if vars.vis == 0, close(figure(160)), end
    
    figure(161)
    imshow3D([f1,f2a,f2b,f3],[0 0.2]), colorbar, title('f'), if vars.saveQA==1, saveas(gcf,fullfile(vars.dir.QApath,'Maps','f'),'fig'), end
    if vars.vis == 0, close(figure(161)), end
    
    figure(162)
    imshow3D([Dstar1,Dstar2a,Dstar2b,Dstar3],[0 .1]), colorbar, title('Dstar'), if vars.saveQA==1, saveas(gcf,fullfile(vars.dir.QApath,'Maps','Dstar'),'fig'), end
    if vars.vis == 0, close(figure(162)), end
    
    figure(163)
    imshow3D([fDstar1,fDstar2a,fDstar2b,fDstar3],[0 .01]), colorbar, title('Dstarf'), if vars.saveQA==1, saveas(gcf,fullfile(vars.dir.QApath,'Maps','Dstarf'),'fig'), end
    if vars.vis == 0, close(figure(163)), end
    
    figure(164), for k = 1:vars.params.slices
        imagesc([D1(:,:,k),D2a(:,:,k),D2b(:,:,k),D3(:,:,k)]), axis image, caxis([0 2.5e-3]), colorbar, title(['D Fit Slice#',num2str(k)]), if vars.saveQA==1, saveas(gcf,fullfile(vars.dir.QApath,'Maps','SliceMaps',['D_slice',num2str(k)]),'fig'), end
        imagesc([f1(:,:,k),f2a(:,:,k),f2b(:,:,k),f3(:,:,k)]), axis image, caxis([0 0.2]), colorbar, title(['f Fit Slice#',num2str(k)]), if vars.saveQA==1, saveas(gcf,fullfile(vars.dir.QApath,'Maps','SliceMaps',['f_slice',num2str(k)]),'fig'), end
        imagesc([Dstar1(:,:,k),Dstar2a(:,:,k),Dstar2b(:,:,k),Dstar3(:,:,k)]), axis image, caxis([0 0.1]), colorbar, title(['Dstar Fit Slice#',num2str(k)]),  if vars.saveQA==1, saveas(gcf,fullfile(vars.dir.QApath,'Maps','SliceMaps',['Dstar_slice',num2str(k)]),'fig'), end
        imagesc([fDstar1(:,:,k),fDstar2a(:,:,k),fDstar2b(:,:,k),fDstar3(:,:,k)]), axis image, caxis([0 0.01]), colorbar, title(['Dstarf Fit Slice#',num2str(k)]),  if vars.saveQA==1, saveas(gcf,fullfile(vars.dir.QApath,'Maps','SliceMaps',['Dstarf_slice',num2str(k)]),'fig'), end
        imagesc([r1(:,:,k),r2a(:,:,k),r2b(:,:,k),r3(:,:,k)]), axis image, caxis([0 0.01]), colorbar, title(['Resnorm Slice#',num2str(k)]),  if vars.saveQA==1, saveas(gcf,fullfile(vars.dir.QApath,'Maps','SliceMaps',['Resnorm_slice',num2str(k)]),'fig'), end
    end, close(figure(164))
    
    
    figure(165), hold on
    histogram(D1(find(mask)),'BinWidth',5e-5)
    histogram(D2a(find(mask)),'BinWidth',5e-5)
    histogram(D2b(find(mask)),'BinWidth',5e-5)
    histogram(D3(find(mask)),'BinWidth',5e-5)
    title('Distribution of D (s/mm^2) in eroded mask')
    legend('D1','D2a','D2b','D3')
    if vars.saveQA >0, saveas(gcf,fullfile(vars.dir.QApath,'Maps','D_histogram'),'fig'), end
    if vars.vis == 0, close(figure(165)), end
    
    figure(166), hold on
    histogram(f1(find(mask)),'BinWidth',0.005)
    histogram(f2a(find(mask)),'BinWidth',0.005)
    histogram(f2b(find(mask)),'BinWidth',0.005)
    histogram(f3(find(mask)),'BinWidth',0.005)
    title('Distribution of f in eroded mask')
    legend('f1','f2a','f2b','f3')
    if vars.saveQA >0, saveas(gcf,fullfile(vars.dir.QApath,'Maps','f_histogram'),'fig'), end
    if vars.vis == 0, close(figure(166)), end
    
    figure(167), hold on
    histogram(Dstar1(find(mask)),'BinWidth',0.005)
    histogram(Dstar2a(find(mask)),'BinWidth',0.005)
    histogram(Dstar2b(find(mask)),'BinWidth',0.005)
    histogram(Dstar3(find(mask)),'BinWidth',0.005)
    title('Distribution of D* (s/mm^2) in eroded mask')
    legend('D*1','D*2a','D*2b','D*3')
    if vars.saveQA>0, saveas(gcf,fullfile(vars.dir.QApath,'Maps','Dstar_histogram'),'fig'), end
    if vars.vis == 0, close(figure(167)), end
    
    figure(168), hold on
    histogram(fDstar1(find(mask)),'BinWidth',0.001)
    histogram(fDstar2a(find(mask)),'BinWidth',0.001)
    histogram(fDstar2b(find(mask)),'BinWidth',0.001)
    histogram(fDstar3(find(mask)),'BinWidth',0.001)
    title('Distribution of D*f (s/mm^2) in eroded mask')
    legend('D*f1','D*f2a','D*f2a','D*f3')
    if vars.saveQA >0, saveas(gcf,fullfile(vars.dir.QApath,'Maps','fDstar_histogram'),'fig'), end
    if vars.vis == 0, close(figure(168)), end
    
end
