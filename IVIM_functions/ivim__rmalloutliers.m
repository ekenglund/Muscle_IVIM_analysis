function ivim__rmalloutliers(iSubj,iPrePost)
disp('Saving voxelwise map data... ')
global vars init_vars
% difference between rmalloutliers_2 and rmalloutliers_1 is that the
% boundary exclusion criteria used for D* is actually the minimum bound,
% whereas in rmalloutliers_1 it is the theoretical boundary (0.00251).
% Also upper bound is much lower in _1 than it should be.. Modified that
% here such that 
%% Extract data from vars


D1 = vars.ivim.All_D(:,:,:,1);
f1 = vars.ivim.All_f(:,:,:,1);
Dstar1 = vars.ivim.All_Dstar(:,:,:,1);
fDstar1 = vars.ivim.All_fDstar(:,:,:,1);
r1 = vars.ivim.All_resnorm(:,:,:,1);

D2a = vars.ivim.All_D(:,:,:,2);
f2a = vars.ivim.All_f(:,:,:,2);
Dstar2a = vars.ivim.All_Dstar(:,:,:,2);
fDstar2a = vars.ivim.All_fDstar(:,:,:,2);
r2a = vars.ivim.All_resnorm(:,:,:,2);


D2b = vars.ivim.All_D(:,:,:,3);
f2b = vars.ivim.All_f(:,:,:,3);
Dstar2b = vars.ivim.All_Dstar(:,:,:,3);
fDstar2b = vars.ivim.All_fDstar(:,:,:,3);
r2b = vars.ivim.All_resnorm(:,:,:,3);


D3 = vars.ivim.All_D(:,:,:,4);
f3 = vars.ivim.All_f(:,:,:,4);
Dstar3 = vars.ivim.All_Dstar(:,:,:,4);
fDstar3 = vars.ivim.All_fDstar(:,:,:,4);
r3 = vars.ivim.All_resnorm(:,:,:,4);

mask = vars.mask.mask_under1;

%% Clean data to ignore points on boundaries
% note that this also exlcudes NaN values
%modified to include boundaries
% D_lim = [-1 1];
% f_lim = [-1 1];
% Dstar_lim = [-1 1];
% fDstar_lim = [-1 1];

%actual limits used
D_lim = [0.00001,0.00249];
f_lim = [0.00001,0.49];
Dstar_lim = [0.00151,0.49];
fDstar_lim = [0.00001,0.1249];

D1_mask = zeros(size(D1));
f1_mask = zeros(size(f1));
Dstar1_mask = zeros(size(Dstar1));
fDstar1_mask = zeros(size(fDstar1));

for r = 1:size(D1,1)
    for c = 1:size(D1,2)
        for s = 1:size(D1,3)
            if D1(r,c,s)>D_lim(1)&&D1(r,c,s)<D_lim(2)
                D1_mask(r,c,s) = 1;
            end
            if f1(r,c,s)>f_lim(1)&&f1(r,c,s)<f_lim(2)
                f1_mask(r,c,s) = 1;
            end
            if Dstar1(r,c,s)>Dstar_lim(1)&&Dstar1(r,c,s)<Dstar_lim(2)
                Dstar1_mask(r,c,s) = 1;
            end
            if fDstar1(r,c,s)>fDstar_lim(1)&&fDstar1(r,c,s)<fDstar_lim(2)
                fDstar1_mask(r,c,s) = 1;
            end
        end
    end
end

for k = 1:size(D1,3)
    mask1(:,:,k) = D1_mask(:,:,k).*f1_mask(:,:,k).*Dstar1_mask(:,:,k).*fDstar1_mask(:,:,k);
end



D2a_mask = zeros(size(D2a));
f2a_mask = zeros(size(f2a));
Dstar2a_mask = zeros(size(Dstar2a));
fDstar2a_mask = zeros(size(fDstar2a));

for r = 1:size(D2a,1)
    for c = 1:size(D2a,2)
        for s = 1:size(D2a,3)
            if D2a(r,c,s)>D_lim(1)&&D2a(r,c,s)<D_lim(2)
                D2a_mask(r,c,s) = 1;
            end
            if f2a(r,c,s)>f_lim(1)&&f2a(r,c,s)<f_lim(2)
                f2a_mask(r,c,s) = 1;
            end
            if Dstar2a(r,c,s)>Dstar_lim(1)&&Dstar2a(r,c,s)<Dstar_lim(2)
                Dstar2a_mask(r,c,s) = 1;
            end
            if fDstar2a(r,c,s)>fDstar_lim(1)&&fDstar2a(r,c,s)<fDstar_lim(2)
                fDstar2a_mask(r,c,s) = 1;
            end
        end
    end
end

for k = 1:size(D2a,3)
    mask2a(:,:,k) = D2a_mask(:,:,k).*f2a_mask(:,:,k).*Dstar2a_mask(:,:,k).*fDstar2a_mask(:,:,k);
end



D2b_mask = zeros(size(D2b));
f2b_mask = zeros(size(f2b));
Dstar2b_mask = zeros(size(Dstar2b));
fDstar2b_mask = zeros(size(fDstar2b));

for r = 1:size(D2b,1)
    for c = 1:size(D2b,2)
        for s = 1:size(D2b,3)
            if D2b(r,c,s)>D_lim(1)&&D2b(r,c,s)<D_lim(2)
                D2b_mask(r,c,s) = 1;
            end
            if f2b(r,c,s)>f_lim(1)&&f2b(r,c,s)<f_lim(2)
                f2b_mask(r,c,s) = 1;
            end
            if Dstar2b(r,c,s)>Dstar_lim(1)&&Dstar2b(r,c,s)<Dstar_lim(2)
                Dstar2b_mask(r,c,s) = 1;
            end
            if fDstar2b(r,c,s)>fDstar_lim(1)&&fDstar2b(r,c,s)<fDstar_lim(2)
                fDstar2b_mask(r,c,s) = 1;
            end
        end
    end
end

for k = 1:size(D2b,3)
    mask2b(:,:,k) = D2b_mask(:,:,k).*f2b_mask(:,:,k).*Dstar2b_mask(:,:,k).*fDstar2b_mask(:,:,k);
end


D3_mask = zeros(size(D3));
f3_mask = zeros(size(f3));
Dstar3_mask = zeros(size(Dstar3));
fDstar3_mask = zeros(size(fDstar3));

for r = 1:size(D3,1)
    for c = 1:size(D3,2)
        for s = 1:size(D3,3)
            if D3(r,c,s)>D_lim(1)&&D3(r,c,s)<D_lim(2)
                D3_mask(r,c,s) = 1;
            end
            if f3(r,c,s)>f_lim(1)&&f3(r,c,s)<f_lim(2)
                f3_mask(r,c,s) = 1;
            end
            if Dstar3(r,c,s)>Dstar_lim(1)&&Dstar3(r,c,s)<Dstar_lim(2)
                Dstar3_mask(r,c,s) = 1;
            end
            if fDstar3(r,c,s)>fDstar_lim(1)&&fDstar3(r,c,s)<fDstar_lim(2)
                fDstar3_mask(r,c,s) = 1;
            end
        end
    end
end

for k = 1:size(D3,3)
    mask3(:,:,k) = D3_mask(:,:,k).*f3_mask(:,:,k).*Dstar3_mask(:,:,k).*fDstar3_mask(:,:,k);
end

vars.mask.fullexclusion.mask1 = mask1;
vars.mask.fullexclusion.mask2a = mask2a;
vars.mask.fullexclusion.mask2b = mask2b;
vars.mask.fullexclusion.mask3 = mask3;
%% 

D1 = D1(find(mask1));
f1 = f1(find(mask1));
Dstar1 = Dstar1(find(mask1));
fDstar1 = fDstar1(find(mask1));
r1 = r1(find(mask1));

D2a = D2a(find(mask2a));
f2a = f2a(find(mask2a));
Dstar2a = Dstar2a(find(mask2a));
fDstar2a = fDstar2a(find(mask2a));
r2a = r2a(find(mask2a));

D2b = D2b(find(mask2b));
f2b = f2b(find(mask2b));
Dstar2b = Dstar2b(find(mask2b));
fDstar2b = fDstar2b(find(mask2b));
r2b = r2b(find(mask2b));

D3 = D3(find(mask3));
f3 = f3(find(mask3));
Dstar3 = Dstar3(find(mask3));
fDstar3 = fDstar3(find(mask3));
r3 = r3(find(mask3));

%%

D1_mean = mean(D1);%(find(mask)));
f1_mean = mean(f1);%(find(mask)));
Dstar1_mean = mean(Dstar1);%(find(mask)));
fDstar1_mean = mean(fDstar1);%(find(mask)));
r1_mean = mean(r1);
fit1(:,1) = [D1_mean; f1_mean; Dstar1_mean; fDstar1_mean; r1_mean];
D1_median = median(D1);%(find(mask)));
f1_median = median(f1);%(find(mask)));
Dstar1_median = median(Dstar1);%(find(mask)));
fDstar1_median = median(fDstar1);%(find(mask)));
r1_median = median(r1);
fit1(:,2) = [D1_median; f1_median; Dstar1_median; fDstar1_median; r1_median];
D1_std = std(D1);%(find(mask)));
f1_std = std(f1);%(find(mask)));
Dstar1_std = std(Dstar1);%(find(mask)));
fDstar1_std = std(fDstar1);%(find(mask)));
r1_std = std(r1);
fit1(:,3) = [D1_std; f1_std; Dstar1_std; fDstar1_std; r1_std];


D2a_mean = mean(D2a);%(find(mask)));
f2a_mean = mean(f2a);%(find(mask)));
Dstar2a_mean = mean(Dstar2a);%(find(mask)));
fDstar2a_mean = mean(fDstar2a);%(find(mask)));
r2a_mean = mean(r2a);
fit2a(:,1) = [D2a_mean; f2a_mean; Dstar2a_mean; fDstar2a_mean; r2a_mean];
D2a_median = median(D2a);%(find(mask)));
f2a_median = median(f2a);%(find(mask)));
Dstar2a_median = median(Dstar2a);%(find(mask)));
fDstar2a_median = median(fDstar2a);%(find(mask)));
r2a_median = median(r2a);
fit2a(:,2) = [D2a_median; f2a_median; Dstar2a_median; fDstar2a_median; r2a_median];
D2a_std = std(D2a);%(find(mask)));
f2a_std = std(f2a);%(find(mask)));
Dstar2a_std = std(Dstar2a);%(find(mask)));
fDstar2a_std = std(fDstar2a);%(find(mask)));
r2a_std = std(r2a);
fit2a(:,3) = [D2a_std; f2a_std; Dstar2a_std; fDstar2a_std; r2a_std];


D2b_mean = mean(D2b);%(find(mask)));
f2b_mean = mean(f2b);%(find(mask)));
Dstar2b_mean = mean(Dstar2b);%(find(mask)));
fDstar2b_mean = mean(fDstar2b);%(find(mask)));
r2b_mean = mean(r2b);
fit2b(:,1) = [D2b_mean; f2b_mean; Dstar2b_mean; fDstar2b_mean; r2b_mean];
D2b_median = median(D2b);%(find(mask)));
f2b_median = median(f2b);%(find(mask)));
Dstar2b_median = median(Dstar2b);%(find(mask)));
fDstar2b_median = median(fDstar2b);%(find(mask)));
r2b_median = median(r2b);
fit2b(:,2) = [D2b_median; f2b_median; Dstar2b_median; fDstar2b_median; r2b_median];
D2b_std = std(D2b);%(find(mask)));
f2b_std = std(f2b);%(find(mask)));
Dstar2b_std = std(Dstar2b);%(find(mask)));
fDstar2b_std = std(fDstar2b);%(find(mask)));
r2b_std = std(r2b);
fit2b(:,3) = [D2b_std; f2b_std; Dstar2b_std; fDstar2b_std; r2b_std];


D3_mean = mean(D3);%(find(mask)));
f3_mean = mean(f3);%(find(mask)));
Dstar3_mean = mean(Dstar3);%(find(mask)));
fDstar3_mean = mean(fDstar3);%(find(mask)));
r3_mean = mean(r3);
fit3(:,1) = [D3_mean; f3_mean; Dstar3_mean; fDstar3_mean; r3_mean];
D3_median = median(D3);%(find(mask)));
f3_median = median(f3);%(find(mask)));
Dstar3_median = median(Dstar3);%(find(mask)));
fDstar3_median = median(fDstar3);%(find(mask)));
r3_median = median(r3);
fit3(:,2) = [D3_median; f3_median; Dstar3_median; fDstar3_median; r3_median];
D3_std = std(D3);%(find(mask)));
f3_std = std(f3);%(find(mask)));
Dstar3_std = std(Dstar3);%(find(mask)));
fDstar3_std = std(fDstar3);%(find(mask)));
r3_std = std(r3);
fit3(:,3) = [D3_std; f3_std; Dstar3_std; fDstar3_std; r3_std];


%         save(fullfile(ResultsPath,'Average_Results_erodedMask_modSlices.mat'),'fit1','fit2','fit3','fit4');


%%
% structure of fit = column 1 --> f, column 2 --> D*, column
% 3 --> D, column 4 --> fD*, column 5 --> resnorm
% row 1 --> method 1, row 2 --> 2a, row 3 --> 2b, row 4 --> 3
vars.ivim.fits_mean = [fit1(2,1), fit1(3,1), fit1(1,1), fit1(4,1), fit1(5,1);...
    fit2a(2,1), fit2a(3,1), fit2a(1,1), fit2a(4,1), fit2a(5,1);...
    fit2b(2,1), fit2b(3,1), fit2b(1,1), fit2b(4,1), fit2b(5,1);...
    fit3(2,1), fit3(3,1), fit3(1,1), fit3(4,1), fit3(5,1)];
% Now structure of this matches the structure of
% ROIAverage_thenFit_Results (plus addition of fD* export)


vars.ivim.fits_median = [fit1(2,2), fit1(3,2), fit1(1,2), fit1(4,2), fit1(5,2);...
    fit2a(2,2), fit2a(3,2), fit2a(1,2), fit2a(4,2), fit2a(5,2);...
    fit2b(2,2), fit2b(3,2), fit2b(1,2), fit2b(4,2), fit2b(5,2);...
    fit3(2,2), fit3(3,2), fit3(1,2), fit3(4,2), fit3(5,2)];


vars.ivim.fits_stdev = [fit1(2,3), fit1(3,3), fit1(1,3), fit1(4,3), fit1(5,3);...
    fit2a(2,3), fit2a(3,3), fit2a(1,3), fit2a(4,3), fit2a(5,3);...
    fit2b(2,3), fit2b(3,3), fit2b(1,3), fit2b(4,3), fit2b(5,3);...
    fit3(2,3), fit3(3,3), fit3(1,3), fit3(4,3), fit3(5,3)];

vars.mask.mask_sizes = [length(find(mask1)); length(find(mask2a)); length(find(mask2b)); length(find(mask3)); length(find(mask))];
% Uncomment below to write results to a file.
% xlwrite(fullfile(vars.dir.basedir,'Summary_Results',['20200709_IVIM_ParameterMapsResults_exmask_rmboundsall_',num2str(vars.MC_option),'.xls']),fits_mean,[init_vars.Subj_list{iSubj},'_',num2str(iPrePost)],'A1:E4');
% xlwrite(fullfile(vars.dir.basedir,'Summary_Results',['20200709_IVIM_ParameterMapsResults_exmask_rmboundsall_',num2str(vars.MC_option),'.xls']),fits_median,[init_vars.Subj_list{iSubj},'_',num2str(iPrePost)],'G1:K4');
% xlwrite(fullfile(vars.dir.basedir,'Summary_Results',['20200709_IVIM_ParameterMapsResults_exmask_rmboundsall_',num2str(vars.MC_option),'.xls']),fits_stdev,[init_vars.Subj_list{iSubj},'_',num2str(iPrePost)],'M1:Q4');
% xlwrite(fullfile(vars.dir.basedir,'Summary_Results',['20200709_IVIM_ParameterMapsResults_exmask_rmboundsall_',num2str(vars.MC_option),'.xls']),mask_sizes,[init_vars.Subj_list{iSubj},'_',num2str(iPrePost)],'S1:S5');
