function [coeffs,resnorms] = ivim__parameter_fitting(SbSo)

global vars

SbSo = double(SbSo);
b_thresh = find(vars.ivim.b>vars.params.highb_threshold,1);
b = vars.ivim.b;

% Fitting will be done 4 ways:
% 1: All parameters are fit from bi-exponential simultaneously
% 2: Fit D from high b-value data (without allowing for an offset at b=0), then fit D* and f from all data, keeping D constant
% 3: Fit D from high b-value data (DO allow for an offset at b=0, but do not use it for anything, then fit D* and f from all data, keeping D constant
% 4: Fit D from log(high b-value data) use intercept to solve for f, then fit for D*, keeping D and f constant


fun_off=inline('(1-ff(1))*exp(-bbbb*ff(2))','ff','bbbb');
funt=inline('f(1)*exp(-bbb(:,1).*f(2))+(1-f(1))*exp(-bbb(:,1).*f(3))','f','bbb');
fun=inline('exp(-bb*dd)','bb','dd');


temp1 = isnan([1 SbSo]);
temp2 = isinf([1 SbSo]);
temp = abs(temp1+temp2-1);
temp_highb = zeros(size(temp));
temp_highb(b_thresh:end) = temp(b_thresh:end);

SbSo_temp = [1 SbSo];
log_SbSo_temp = log(SbSo_temp);


if sum(temp)>3 % there must be at least three non-zero values to fit with 1 part fitting
    % fit method 1
    [coeffs_1partfit,resnorm_1partfit,residual_1partfit] = ...
        lsqcurvefit(funt,[0.1 0.01 0.0015],b(find(temp))',SbSo_temp(find(temp))',[0 .0015 0],...
        [.5 .5 0.0025],optimset('Display','off'));
else
    coeffs_1partfit=[NaN,NaN,NaN];
    resnorm_1partfit = NaN;
end

% fit methods 2-4 only work if there is enough high b-value data
if sum(temp_highb)>2
    % fit method 2
    coeff_D_highb_0=lsqcurvefit(fun,0.0015,b(find(temp_highb)),...
        SbSo_temp(find(temp_highb)),0,0.0025,optimset('Display','off'));
    [coeffs_2partfita,resnorm_2partfita,residual_2partfita] = ...
        lsqcurvefit(funt,[0.1 0.01 coeff_D_highb_0],b(find(temp))',SbSo_temp(find(temp))',...
        [0 .0015 coeff_D_highb_0-1e-6],[.5 .5 coeff_D_highb_0+1e-6],optimset('Display','off'));
    
    % fit method 3
    coeff_D_highb_f=lsqcurvefit(fun_off,[0 0.0015],b(find(temp_highb)),...
        SbSo_temp(find(temp_highb)),[0 0],[0.5 0.0025],optimset('Display','off'));
    [coeffs_2partfitb,resnorm_2partfitb,residual_2partfitb] = ...
        lsqcurvefit(funt,[coeff_D_highb_f(1) 0.01 coeff_D_highb_f(2)],b(find(temp))',...
        SbSo_temp(find(temp))',[0 .0015 coeff_D_highb_f(2)-1e-6],...
        [.5 .5 coeff_D_highb_f(2)+1e-6],optimset('Display','off'));
    
    % fit method 4
    % In full Paraspinal ROI
    Fit= polyfit(b(find(temp_highb)),log_SbSo_temp(find(temp_highb)),1);
    SoPrime = exp(Fit(2));
    D_3partfit = -Fit(1);
    f_3partfit = 1-SoPrime;
    if D_3partfit >2.5e-3
        D_3partfit = 2.5e-3;
    elseif D_3partfit <0
        D_3partfit = 0;
    end
    f_3partfit = 1-SoPrime;
    if f_3partfit < 0
        f_3partfit = 0;
    elseif f_3partfit > 0.5
        f_3partfit = 0.5;
    end
    [coeffs_3partfit,resnorm_3partfit,residual_3partfit] = lsqcurvefit(funt,...
        [f_3partfit 0.01 D_3partfit],b(find(temp))',SbSo_temp(find(temp))',...
        [f_3partfit-1e-6 .0015 D_3partfit-1e-6],[f_3partfit+1e-6 .5 D_3partfit+1e-6],optimset('Display','off'));
else
    coeffs_2partfita = [NaN,NaN,NaN];
    coeffs_2partfitb = [NaN,NaN,NaN];
    coeffs_3partfit = [NaN,NaN,NaN];
    resnorm_2partfita = NaN;
    resnorm_2partfitb = NaN;
    resnorm_3partfit = NaN;
end
coeffs = [coeffs_1partfit;coeffs_2partfita;coeffs_2partfitb;coeffs_3partfit];
resnorms = [resnorm_1partfit; resnorm_2partfita; resnorm_2partfitb; resnorm_3partfit];
% residuals = [residual_1partfit; residual_2partfita; residual_2partfitb; residual_3partfit];



