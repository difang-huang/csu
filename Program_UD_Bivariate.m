%%%%%%%%%%%%%%%%
% Program_main.m
%%%%%%%%%%%%%%%%

% Last modified: 23-05-2021

clear;clc;


% Include stuff for writing to Excel file (if using Mac)

javaaddpath('poi_library/poi-3.8-20120326.jar');
javaaddpath('poi_library/poi-ooxml-3.8-20120326.jar');
javaaddpath('poi_library/poi-ooxml-schemas-3.8-20120326.jar');
javaaddpath('poi_library/xmlbeans-2.3.0.jar');
javaaddpath('poi_library/dom4j-1.6.1.jar');
javaaddpath('poi_library/stax-api-1.0.1.jar');

% Load equity risk premium data, 1996:01-2019:08

input_file='PredictorData2019.xlsx';
input_sheet='Monthly';
Rfree_lag=xlsread(input_file,input_sheet,'k1502:k1785');
R_SP500=xlsread(input_file,input_sheet,'q1502:q1785');
r=log(1+R_SP500)-log(1+Rfree_lag);
disp('Equity risk premium, ann % summary stats (mean, vol, Sharpe ratio)');
disp([12*mean(r) sqrt(12)*std(r) sqrt(12)*mean(r)/std(r)]);

% Load Goyal-Welch predictor data and csu data, 1996:01-2019:08
input_file='csu_monthly.xlsx';
input_sheet='csu_monthly';
EWSI=xlsread(input_file,input_sheet,'b194:b477');
log_EWSI = log(EWSI);
SII = zscore(EWSI);

load('Program_generate_UD_predictors.mat');
stats_GW=[mean(GW)' median(GW)' prctile(GW,1)' prctile(GW,99)' std(GW)'];
disp('Predictor variables, summary stats');
disp('Mean, median, 1st percentile, 99th percentile, std dev');
disp(stats_GW);
rho_GW=nan(size(GW,2));
for i=1:size(GW,2);
    rho_GW(i)=corr(GW(2:end,i),GW(1:end-1,i));
end;
GW_adjust=GW;
GW_adjust(:,7:9)=-GW(:,7:9);
GW_adjust(:,end)=-GW(:,end);
GW_standardize=zscore(GW_adjust);

% Compute principal components for Goyal-Welch predictors

% X=GW_standardize;
% X(:,[4 11])=[];
% [coeff,score,latent]=pca(X);
% PC_GW=zscore(score(:,1:3));

 

% Compute cumulative returns

h=[1 3 6 12];
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;

r_h(isnan(r_h))=0;



% Take care of out-of-sample preliminaries

T = length(r);
in_sample_end = 2000;
R = (in_sample_end-1996)*12; % in-sample period
P = T-R; % out-of-sample period
FC_PM = nan(P,1);
FC_PR = nan(P,size(GW,2),length(h));
FC_SII = nan(P,1,length(h));
FC_OTHER1 = nan(P,size(GW,2),length(h));
FC_OTHER2 = nan(P,size(GW,2),length(h));

% Compute out-of-sample forecasts

for p=1:P;
    disp(p);

    % Prevailing mean benchmark forecast

    FC_PM(p)=mean(r(1:R+(p-1)));

    % Predictive regression forecasts

    for j=1:length(h);
        
        % SII 
        X_SII_j_p=[ones(R+(p-1)-h(j),1) SII(1:R+(p-1)-h(j))];
        results_SII_j_p=ols(r_h(2:R+p-h(j),j),X_SII_j_p);
        FC_SII(p,1,j)=[1 SII(R+(p-1))]*results_SII_j_p.beta;

        % Goyal-Welch predictors

        for i=1:size(GW,2);
            X_i_j_p=[ones(R+(p-1)-h(j),1) GW(1:R+(p-1)-h(j),i)];
            results_i_j_p=ols(r_h(2:R+p-h(j),j),X_i_j_p);
            FC_PR(p,i,j)=[1 GW(R+(p-1),i)]*results_i_j_p.beta;
            FC_OTHER1(p,i,j)= 0.5*FC_PR(p,i,j) + 0.5*FC_SII(p,1,j);
            FC_OTHER2(p,i,j)= 0.2*FC_PR(p,i,j) + 0.8*FC_SII(p,1,j);
        end;

    end;
end;



% Evaluate forecasts

R2OS_PR=nan(size(GW,2),2,length(h));
R2OS_OTHER1=nan(size(GW,2),2,length(h));
R2OS_OTHER2=nan(size(GW,2),2,length(h));

for j=1:length(h);
    actual_j=r_h(R+1:end-(h(j)-1),j);
    u_PM_j=actual_j-FC_PM(1:end-(h(j)-1));
    u_PR_j=kron(ones(1,size(FC_PR,2)),actual_j)-FC_PR(1:end-(h(j)-1),:,j);
    u_PR1_j=kron(ones(1,size(FC_OTHER1,2)),actual_j)-FC_OTHER1(1:end-(h(j)-1),:,j);
    u_PR2_j=kron(ones(1,size(FC_OTHER2,2)),actual_j)-FC_OTHER2(1:end-(h(j)-1),:,j);
    MSFE_PM_j=mean(u_PM_j.^2);
    MSFE_PR_j=mean(u_PR_j.^2);
    MSFE_PR1_j=mean(u_PR1_j.^2);
    MSFE_PR2_j=mean(u_PR2_j.^2);
    R2OS_PR_j=100*(1-MSFE_PR_j/MSFE_PM_j);
    R2OS_PR1_j=100*(1-MSFE_PR1_j/MSFE_PM_j);
    R2OS_PR2_j=100*(1-MSFE_PR2_j/MSFE_PM_j);
    R2OS_PR(:,1,j)=R2OS_PR_j';
    R2OS_OTHER1(:,1,j)=R2OS_PR1_j';
    R2OS_OTHER2(:,1,j)=R2OS_PR2_j';
    L_j=h(j);
    for i=1:size(GW,2);
        f_CW_i_j=u_PM_j.^2-u_PR_j(:,i).^2+(FC_PM(1:end-(h(j)-1))-...
            FC_PR(1:end-(h(j)-1),i,j)).^2;
        results_CW_i_j=nwest(f_CW_i_j,ones(length(f_CW_i_j),1),h(j));
        R2OS_PR(i,2,j) = (1-tcdf(abs(results_CW_i_j.tstat),results_CW_i_j.nobs-2));
        
        f_CW_i_j=u_PM_j.^2-u_PR1_j(:,i).^2+(FC_PM(1:end-(h(j)-1))-...
            FC_OTHER1(1:end-(h(j)-1),i,j)).^2;
        results_CW_i_j=nwest(f_CW_i_j,ones(length(f_CW_i_j),1),h(j));
        R2OS_OTHER1(i,2,j)= (1-tcdf(abs(results_CW_i_j.tstat),results_CW_i_j.nobs-2));
        
        f_CW_i_j=u_PM_j.^2-u_PR2_j(:,i).^2+(FC_PM(1:end-(h(j)-1))-...
            FC_OTHER2(1:end-(h(j)-1),i,j)).^2;
        results_CW_i_j=nwest(f_CW_i_j,ones(length(f_CW_i_j),1),h(j));
        R2OS_OTHER2(i,2,j)= (1-tcdf(abs(results_CW_i_j.tstat),results_CW_i_j.nobs-2));
    end;
end;

output_file='Results_UD.xlsx';

output_sheet='R2OS_OTHER1 statistics';
xlwrite(output_file,R2OS_OTHER1(:,:,1),output_sheet,'b3');
xlwrite(output_file,R2OS_OTHER1(:,:,2),output_sheet,'e3');
xlwrite(output_file,R2OS_OTHER1(:,:,3),output_sheet,'h3');
xlwrite(output_file,R2OS_OTHER1(:,:,4),output_sheet,'k3');

output_sheet='R2OS_OTHER2 statistics';
xlwrite(output_file,R2OS_OTHER2(:,:,1),output_sheet,'b3');
xlwrite(output_file,R2OS_OTHER2(:,:,2),output_sheet,'e3');
xlwrite(output_file,R2OS_OTHER2(:,:,3),output_sheet,'h3');
xlwrite(output_file,R2OS_OTHER2(:,:,4),output_sheet,'k3');
 
