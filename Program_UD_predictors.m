%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program_generate_GW_predictors.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Last modified: 08-18-2015

clear;

% Load uncertainty data, 1996:01-2019:08
input_file='UncertaintyData2019.xlsx';
input_sheet='Month';

% (1) VRP

VRP=xlsread(input_file,input_sheet,'c74:c357');

% (2) EVRP

EVRP=xlsread(input_file,input_sheet,'d74:d357');

% (3) IV

IV=xlsread(input_file,input_sheet,'e74:e357');

% (4) RV

RV=xlsread(input_file,input_sheet,'f74:f357');

% (5) ERV

ERV=xlsread(input_file,input_sheet,'g74:g357');

% (6) EUI

EUI=xlsread(input_file,input_sheet,'h74:h357');

% (7) Liquidity

LIQ=xlsread(input_file,input_sheet,'i74:i357');

% (8) CATFIN

CATFIN=xlsread(input_file,input_sheet,'j74:j357');

% (9) BBD

BBD=xlsread(input_file,input_sheet,'k74:k357');

% (10) JLN

JLN=xlsread(input_file,input_sheet,'l74:l357');


% Collect economics variables, 1996:01-2019:08

GW=[VRP EVRP IV RV ERV EUI LIQ CATFIN BBD JLN];
save('Program_generate_UD_predictors.mat','GW');
