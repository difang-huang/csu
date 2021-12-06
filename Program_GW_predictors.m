%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program_generate_GW_predictors.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Last modified: 08-18-2015

clear;

% Load Goyal-Welch predictor data, 1996:01-2019:08
input_file='PredictorData2019.xlsx';
input_sheet='Monthly';

% (1) Log dividend-price ratio

SP=xlsread(input_file,input_sheet,'b1502:b1785');
D12=xlsread(input_file,input_sheet,'c1502:c1785');
log_DP=log(D12./SP);

% (2) Log dividend yield

SP_lag=xlsread(input_file,input_sheet,'b1501:b1784');
log_DY=log(D12./SP_lag);

% (3) Log earnings-price ratio

E12=xlsread(input_file,input_sheet,'d1502:d1785');
log_EP=log(E12./SP);

% (4) Log dividend-payout ratio

log_DE=log(D12./E12);

% (5) Stock excess return volatility (annualized)

SP_R=xlsread(input_file,input_sheet,'q1491:q1785'); % 1995:02-2019:08
R_F_lag=xlsread(input_file,input_sheet,'k1490:k1784'); % 1995:01-2019:07
RVOL=nan(length(SP_R)-11,1); % using Mele (2007) estimator
for t=1:length(RVOL);
    RVOL(t)=mean(abs(SP_R(t:t+11)-R_F_lag(t:t+11)));
end;
RVOL=sqrt(pi/2)*sqrt(12)*RVOL; % 1973:01-2014:12

% (6) Book-to-market ratio

BM=xlsread(input_file,input_sheet,'e1502:e1785');

% (7) Net equity issuance

NTIS=xlsread(input_file,input_sheet,'j1502:j1785');

% (8) Treasury bill rate (annualized)

TBL=xlsread(input_file,input_sheet,'f1502:f1785');

% (9) Long-term yield (annualized)

LTY=xlsread(input_file,input_sheet,'i1502:i1785');

% (10) Long-term return

LTR=xlsread(input_file,input_sheet,'m1502:m1785');

% (11) Term spread (annualized)

TMS=LTY-TBL;

% (12) Default yield spread (annualized)

BAA=xlsread(input_file,input_sheet,'h1502:h1785');
AAA=xlsread(input_file,input_sheet,'g1502:g1785');
DFY=BAA-AAA;

% (13) Default return spread

CORPR=xlsread(input_file,input_sheet,'n1502:n1785');
DFR=CORPR-LTR;

% (14) Inflation (lagged, to account for data release)

INFL_lag=xlsread(input_file,input_sheet,'l1501:l1784');

% Collect economics variables, 1973:01-2014:12

GW=[log_DP log_DY log_EP log_DE RVOL BM NTIS TBL LTY LTR TMS DFY DFR ...
    INFL_lag];
save('Program_generate_GW_predictors.mat','GW');
