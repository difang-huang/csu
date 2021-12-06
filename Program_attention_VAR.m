% This program examines the time series data on the csuindex and the
% funding risk measure of Fontaine and Garcia (2012).

clear;clc;

%--------------------------------------------------------------------------
%   LOAD THE DATA
%--------------------------------------------------------------------------

load csu_and_attention;  % Monthly levels of the two indices from
                                  % Dec 1985 to Dec 2011 (T=313)

FR = FR-mean(FR);
Tail = Tail-mean(Tail);

%--------------------------------------------------------------------------
%   SOME PRELIMINARY CALCULATIONS
%--------------------------------------------------------------------------

formatSpec = '%9.3f';
disp('Summary Statistics on the Levels of the Series:');
disp(' ');
disp(['Name     ' 'Mean    ' 'StdDev    ' 'AR(1)    ' 'AR(2)    ' 'AR(3)    ' 'AR(4)    ' 'AR(5)    ' 'AR(6)']); 
disp(['----     ' '----    ' '------    ' '-----    ' '-----    ' '-----    ' '-----    ' '-----    ' '-----']);
disp(['Tail     ' num2str([mean(Tail)  std(Tail) corr(Tail(2:end),Tail(1:end-1)) corr(Tail(3:end),Tail(1:end-2)) ...
    corr(Tail(4:end),Tail(1:end-3)) corr(Tail(5:end),Tail(1:end-4)) corr(Tail(6:end),Tail(1:end-5)) ...
    corr(Tail(7:end),Tail(1:end-6))],formatSpec)]);
disp(['FundMeas ' num2str([mean(FR)  std(FR) corr(FR(2:end),FR(1:end-1)) corr(FR(3:end),FR(1:end-2)) ...
    corr(FR(4:end),FR(1:end-3)) corr(FR(5:end),FR(1:end-4)) corr(FR(6:end),FR(1:end-5)) ...
    corr(FR(7:end),FR(1:end-6))],formatSpec)]);
disp(' ');
disp(' ');

% Univariate tests for I(1) process ...

[adfstat_tail,pval_tail,critval_tail,resid_tail,lags_tail]=augdfautolag(Tail,1,24);
[adfstat_FR,pval_FR,critval_FR,resid_FR,lags_FR]=augdfautolag(FR,1,24);

disp('Recall, H0 of ADF is that the series HAS a unit root ...');
disp(['ADF Test for TAIL = ' num2str(adfstat_tail,formatSpec) '  P-Val = ' num2str(pval_tail,formatSpec)  '  Lags = ' num2str(lags_tail,formatSpec)]);
disp(['ADF Test for FR  = ' num2str(adfstat_FR,formatSpec)  '  P-Val = ' num2str(pval_FR,formatSpec)   '  Lags = ' num2str(lags_FR,formatSpec)]);
disp(' ');

y = [Tail FR];

%--------------------------------------------------------------------------
%   ESTIMATING THE VAR (and Granger Causality): 1st Order VAR
%--------------------------------------------------------------------------

% Estimating the univariate autoregressions
[paramT1,stderrT1,tstatT1,pvalT1,constT1,conststdT1,r2_T1,errorsT1,s2_T1,paramvecT1,vcvT1]  = vectorar(Tail,0,1);
[paramFR1,stderrFR1,tstatFR1,pvalFR1,constFR1,conststdFR1,r2_FR1,errorsFR1,s2_FR1,paramvecFR1,vcvFR1]  = vectorar(FR,0,1);

% Estimating the VAR
[param1,stderr1,tstat1,pval1,const1,conststd1,r2_1,errors1,s2_1,paramvec1,vcv1] = vectorar(y,0,1);

disp('*************************');
disp('First-Order VAR in Levels');
disp('-------------------------');
disp('Coefficients:');
disp(num2str(cell2mat(param1)));
disp('Asympt. t-stats:');
disp(num2str(cell2mat(tstat1)));
disp('');


a = cell2mat(param1);
b = cell2mat(tstat1);
output_file='Results_CSU.xlsx';
output_sheet='Attention VAR';
xlwrite(output_file,a,output_sheet,'b3');
xlwrite(output_file,b,output_sheet,'e3');

formatSpec = '%10.4f';
% Granger Causality Tests in Levels -- implemented following the approach
% in Section 11.2 of Hamilton (specifically, eq. 11.2.10)

%---------------------------------------------
% Does Tail Granger cause FR Index?
RSS1 = errors1(:,2)'*errors1(:,2);
RSS0 = errorsFR1'*errorsFR1;
Test11   = (length(errors1)*(RSS0-RSS1))/RSS1';
Test11pv = 1-chi2cdf(Test11,1);
disp(['First-Order Granger Causality of FR Index by Tail = ' num2str(Test11,formatSpec) '  Asymp. P-Val = ' num2str(Test11pv,formatSpec)]);

%---------------------------------------------
% Does FR Index Granger cause Tail?
RSS1 = errors1(:,1)'*errors1(:,1);
RSS0 = errorsT1'*errorsT1;
Test12   = (length(errors1)*(RSS0-RSS1))/RSS1';
Test12pv = 1-chi2cdf(Test12,1);
disp(['First-Order Granger Causality of Tail by FR Index = ' num2str(Test12,formatSpec) '  Asymp. P-Val = ' num2str(Test12pv,formatSpec)]);
disp(' ');
disp(' ');

    %--------------------------------------------------------
    % Bootstrapping p-values for VAR coefficients (IN LEVELS)
    %--------------------------------------------------------

    % ********* BOOTSTRAPPING THE 1-MONTH VAR *********
    
    param1      = cell2mat(param1);
    BootNum     = 20000;
    BootMat1    = zeros(BootNum,4);
    BootGC1test = zeros(BootNum,2); % Col 1: TAIL ==> FR
                                    % Col 2: FR  ==> TAIL
    [T1,~]      = size(errors1);

    [~, BootIndices1]=stationary_bootstrap(errors1(:,1),BootNum,T1/20);
    
    ii = 1000;
    for i=1:BootNum
    
        zzz1   = mean(errors1(BootIndices1(:,i),:));
        tempE1 = errors1(BootIndices1(:,i),:) - ...
                 ones(size(errors1(BootIndices1(:,i),:)))*diag(zzz1);
 
        % Building the simulated VAR series in each bootstrap sample. Note
        % that I am building each series under the null of no Granger
        % causality. So, I am actually going to generate each series
        % separately.
        
        yboot = zeros(size(y));
        for j=2:T1
            yboot(j,1) = param1(1,1)*yboot(j-1,1)+ tempE1(j,1);
            yboot(j,2) = param1(2,2)*yboot(j-1,2)+ tempE1(j,2);
        end     
        
       [~,~,boot_tstat1,~,~,~,~,booterr1,~,~,~] = vectorar(yboot,0,1);
       boot_t1 = vec(cell2mat(boot_tstat1));
       BootMat1(i,:) = boot_t1';
       [~,~,~,~,~,~,~,booterrFR1,~,~,~]  = vectorar(yboot(:,2),0,1);    
       [~,~,~,~,~,~,~,booterrT1,~,~,~]   = vectorar(yboot(:,1),0,1);
       
       bootRSS1 = booterr1(:,2)'*booterr1(:,2);
       bootRSS0 = booterrFR1'*booterrFR1;
       BootGC1test(i,1) = (length(booterr1)*(bootRSS0-bootRSS1))/bootRSS1';
       
       bootRSS1 = booterr1(:,1)'*booterr1(:,1);
       bootRSS0 = booterrT1'*booterrT1;
       BootGC1test(i,2) = (length(booterr1)*(bootRSS0-bootRSS1))/bootRSS1';
       
       if i == ii
         disp(['Just finished Bootstrap for Levels VAR Length=1:  ' int2str(i) ' of ' int2str(BootNum)]);
         ii = ii+1000;
       end
    end
    
% Computing p-Values (under the null) from the bootstrap ...

p1boot   = sort(abs(BootMat1(:,1)));
p2boot   = sort(abs(BootMat1(:,2)));
p3boot   = sort(abs(BootMat1(:,3)));
p4boot   = sort(abs(BootMat1(:,4)));
pGC11boot = sort(BootGC1test(:,1));
pGC12boot = sort(BootGC1test(:,2));

tstat1 = vec(cell2mat(tstat1));

t1temp = 1-min(find(p1boot>=abs(tstat1(1))))/BootNum;
if isempty(t1temp)==1
    t1temp = 0;
end
t2temp = 1-min(find(p2boot>=abs(tstat1(2))))/BootNum;
if isempty(t2temp)==1
    t2temp = 0;
end
t3temp = 1-min(find(p3boot>=abs(tstat1(3))))/BootNum;
if isempty(t3temp)==1
    t3temp = 0;
end
t4temp = 1-min(find(p4boot>=abs(tstat1(4))))/BootNum;
if isempty(t4temp)==1
    t4temp = 0;
end

GC11temp = 1-min(find(pGC11boot>=Test11))/BootNum;
GC12temp = 1-min(find(pGC12boot>=Test12))/BootNum;

disp(' ');
disp('BOOTSTRAP P-VALUES FROM VAR WITH 1 LAG:');
disp(['VAR coefficients:  ' num2str(vec(param1)')]);
disp(['Bootstrap p-vals:  ' num2str([t1temp t2temp t3temp t4temp])]);
disp(['First-Order Granger Causality of FR Index by Tail = ' num2str(Test11,formatSpec) '  Boot P-Val = ' num2str(GC11temp,formatSpec)]);
disp(['First-Order Granger Causality of Tail Index by FR = ' num2str(Test12,formatSpec) '  Boot P-Val = ' num2str(GC12temp,formatSpec)]);
disp(' ');
disp(' ');

%--------------------------------------------------------------------------
%   ESTIMATING THE VAR (and Granger Causality): 2nd Order VAR 
%--------------------------------------------------------------------------

% Estimating the univariate autoregressions
[paramT2,stderrT2,tstatT2,pvalT2,constT2,conststdT2,r2_T2,errorsT2,s2_T2,paramvecT2,vcvT2]  = vectorar(Tail,0,(1:2));
[paramFR2,stderrFR2,tstatFR2,pvalFR2,constFR2,conststdFR2,r2_FR2,errorsFR2,s2_FR2,paramvecFR2,vcvFR2]  = vectorar(FR,0,(1:2));
% Estimating the VAR
[param2,stderr2,tstat2,pval2,const2,conststd2,r2_2,errors2,s2_2,paramvec2,vcv2] = vectorar(y,0,(1:2));

%---------------------------------------------
% Does Tail Granger cause FR Index?
RSS1 = errors2(:,2)'*errors2(:,2);
RSS0 = errorsFR2'*errorsFR2;
Test21   = (length(errors2)*(RSS0-RSS1))/RSS1';
Test21pv = 1-chi2cdf(Test21,2);
disp(['Second-Order Granger Causality of FR Index by Tail = ' ...
      num2str(Test21,formatSpec) '  Asymp. P-Val = ' num2str(Test21pv,formatSpec)]);

%---------------------------------------------
% Does FR Index Granger cause Tail?
RSS1 = errors2(:,1)'*errors2(:,1);
RSS0 = errorsT2'*errorsT2;
Test22  = (length(errors2)*(RSS0-RSS1))/RSS1';
Test22pv = 1-chi2cdf(Test22,2);
disp(['Second-Order Granger Causality of Tail by FR Index = ' ...
       num2str(Test22,formatSpec) '  Asymp. P-Val = ' num2str(Test22pv,formatSpec)]);
disp(' ');
disp(' ');


    %=====  BOOTSTRAPPING THE SECOND ORDER GRANGER CAUSALITY TEST =====
        
    param2      = cell2mat(param2);
    BootGC2test = zeros(BootNum,2);
    [T2,~]      = size(errors2);

    [~, BootIndices2]=stationary_bootstrap(errors2(:,1),BootNum,T2/20);
    
    ii = 1000;
    for i=1:BootNum
    
        zzz2   = mean(errors2(BootIndices2(:,i),:));
        tempE2 = errors2(BootIndices2(:,i),:) - ...
                 ones(size(errors2(BootIndices2(:,i),:)))*diag(zzz2);
 
        % Building the simulated VAR series in each bootstrap sample. Note
        % that I am building each series under the null of no Granger
        % causality. So, I am actually going to generate each series
        % separately.
        
        yboot2 = zeros(size(y));
        for j=3:T2
            yboot2(j,1) = param2(1,1)*yboot2(j-1,1) + param2(3,1)*yboot2(j-2,1) + tempE2(j,1);
            yboot2(j,2) = param2(2,2)*yboot2(j-1,2) + param2(4,2)*yboot2(j-2,2) + tempE2(j,2);
        end     

       [~,~,~,~,~,~,~,booterrT2,~,~,~]  = vectorar(yboot2(:,1),0,(1:2));
       [~,~,~,~,~,~,~,booterrFR2,~,~,~] = vectorar(yboot2(:,2),0,(1:2));
       [~,~,~,~,~,~,~,booterr2,~,~,~]   = vectorar(yboot2,0,(1:2));
       
       bootRSS1 = booterr2(:,2)'*booterr2(:,2);
       bootRSS0 = booterrFR2'*booterrFR2;
       BootGC2test(i,1) = (length(booterr2)*(bootRSS0-bootRSS1))/bootRSS1';
       bootRSS1 = booterr2(:,1)'*booterr2(:,1);
       bootRSS0 = booterrT2'*booterrT2;
       BootGC2test(i,2) = (length(booterr2)*(bootRSS0-bootRSS1))/bootRSS1';       
       
       if i == ii
         disp(['Just finished Bootstrap for Levels VAR Length=2:  ' int2str(i) ' of ' int2str(BootNum)]);
         ii = ii+1000;
       end
    end

    pGC21boot = sort(BootGC2test(:,1));
    pGC22boot = sort(BootGC2test(:,2));
    
    GC21temp = 1-min(find(pGC21boot>=Test21))/BootNum;
    GC22temp = 1-min(find(pGC22boot>=Test22))/BootNum;
    
    disp(['Second-Order Granger Causality of FR Index by Tail = ' num2str(Test21,formatSpec) '  Boot P-Val = ' num2str(GC21temp,formatSpec)]);
    disp(['Second-Order Granger Causality of Tail Index by FR = ' num2str(Test22,formatSpec) '  Boot P-Val = ' num2str(GC22temp,formatSpec)]);
    disp(' ');
    disp(' ');
    
%--------------------------------------------------------------------------
%    ESTIMATING THE VAR (and Granger Causality): 3rd Order VAR 
%--------------------------------------------------------------------------

% Estimating the univariate autoregressions
[paramT3,stderrT3,tstatT3,pvalT3,constT3,conststdT3,r2_T3,errorsT3,s2_T3,paramvecT3,vcvT3]  = vectorar(Tail,0,(1:3));
[paramFR3,stderrFR3,tstatFR3,pvalFR3,constFR3,conststdFR3,r2_FR3,errorsFR3,s2_FR3,paramvecFR3,vcvFR3]  = vectorar(FR,0,(1:3));
% Estimating the VAR
[param3,stderr3,tstat3,pval3,const3,conststd3,r2_3,errors3,s2_3,paramvec3,vcv3] = vectorar(y,0,(1:3));

%---------------------------------------------
% Does Tail Granger cause FR Index?
RSS1 = errors3(:,2)'*errors3(:,2);
RSS0 = errorsFR3'*errorsFR3;
Test31   = (length(errors3)*(RSS0-RSS1))/RSS1';
Test31pv = 1-chi2cdf(Test31,3);
disp(['Third-Order Granger Causality of FR Index by Tail  = ' num2str(Test31,formatSpec) '  Asymp. P-Val = ' num2str(Test31pv,formatSpec)]);

%---------------------------------------------
% Does FR Index Granger cause Tail?
RSS1 = errors3(:,1)'*errors3(:,1);
RSS0 = errorsT3'*errorsT3;
Test32   = (length(errors3)*(RSS0-RSS1))/RSS1';
Test32pv = 1-chi2cdf(Test32,3);
disp(['Third-Order Granger Causality of Tail by FR Index  = ' num2str(Test32,formatSpec) '  Asymp. P-Val = ' num2str(Test32pv,formatSpec)]);
disp(' ');
disp(' ');


    %=====  BOOTSTRAPPING THE THIRD ORDER GRANGER CAUSALITY TEST =====
        
    param3      = cell2mat(param3);
    BootGC3test = zeros(BootNum,2);
    [T3,~]      = size(errors3);

    [~, BootIndices3]=stationary_bootstrap(errors3(:,1),BootNum,T3/20);
    
    ii = 1000;
    for i=1:BootNum
    
        zzz3   = mean(errors3(BootIndices3(:,i),:));
        tempE3 = errors3(BootIndices3(:,i),:) - ...
                 ones(size(errors3(BootIndices3(:,i),:)))*diag(zzz3);
 
        % Building the simulated VAR series in each bootstrap sample. Note
        % that I am building each series under the null of no Granger
        % causality. So, I am actually going to generate each series
        % separately.
        
        yboot3 = zeros(size(y));
        for j=4:T3
            yboot3(j,1) = param3(1,1)*yboot3(j-1,1) + ...
                          param3(3,1)*yboot3(j-2,1) + ...
                          param3(5,1)*yboot3(j-2,1) + tempE3(j,1);
            yboot3(j,2) = param3(2,2)*yboot3(j-1,2) + ...
                          param3(4,2)*yboot3(j-2,2) + ...
                          param3(6,2)*yboot3(j-2,2) + tempE3(j,2);
        end     

       [~,~,~,~,~,~,~,booterrT3,~,~,~]  = vectorar(yboot3(:,1),0,(1:3));
       [~,~,~,~,~,~,~,booterrFR3,~,~,~] = vectorar(yboot3(:,2),0,(1:3));
       [~,~,~,~,~,~,~,booterr3,~,~,~]   = vectorar(yboot3,0,(1:3));
       
       bootRSS1 = booterr3(:,2)'*booterr3(:,2);
       bootRSS0 = booterrFR3'*booterrFR3;
       BootGC3test(i,1) = (length(booterr3)*(bootRSS0-bootRSS1))/bootRSS1';
       bootRSS1 = booterr3(:,1)'*booterr3(:,1);
       bootRSS0 = booterrT3'*booterrT3;
       BootGC3test(i,2) = (length(booterr3)*(bootRSS0-bootRSS1))/bootRSS1';       
       
       if i == ii
         disp(['Just finished Bootstrap for Levels VAR Length=3  ' int2str(i) ' of ' int2str(BootNum)]);
         ii = ii+1000;
       end
    end

    pGC31boot = sort(BootGC3test(:,1));
    pGC32boot = sort(BootGC3test(:,2));
    
    GC31temp = 1-min(find(pGC31boot>=Test31))/BootNum;
    GC32temp = 1-min(find(pGC32boot>=Test32))/BootNum;
    
    disp(['Third-Order Granger Causality of FR Index by Tail = ' num2str(Test31,formatSpec) '  Boot P-Val = ' num2str(GC31temp,formatSpec)]);
    disp(['Third-Order Granger Causality of Tail Index by FR = ' num2str(Test32,formatSpec) '  Boot P-Val = ' num2str(GC32temp,formatSpec)]);
    disp(' ');
    disp(' ');

    
    
output_file='Results_CSU.xlsx';
output_sheet='Attention VAR';
xlwrite(output_file,Test11,output_sheet,'b6');
xlwrite(output_file,Test12,output_sheet,'b10');
xlwrite(output_file,GC11temp,output_sheet,'b7');
xlwrite(output_file,GC12temp,output_sheet,'b11');
xlwrite(output_file,Test11pv,output_sheet,'b8');
xlwrite(output_file,Test12pv,output_sheet,'b12');


xlwrite(output_file,Test21,output_sheet,'d6');
xlwrite(output_file,Test22,output_sheet,'d10');
xlwrite(output_file,GC21temp,output_sheet,'d7');
xlwrite(output_file,GC22temp,output_sheet,'d11');
xlwrite(output_file,Test21pv,output_sheet,'d8');
xlwrite(output_file,Test22pv,output_sheet,'d12');


xlwrite(output_file,Test31,output_sheet,'g6');
xlwrite(output_file,Test32,output_sheet,'g10');
xlwrite(output_file,GC31temp,output_sheet,'g7');
xlwrite(output_file,GC32temp,output_sheet,'g11');
xlwrite(output_file,Test31pv,output_sheet,'g8');
xlwrite(output_file,Test32pv,output_sheet,'g12');
