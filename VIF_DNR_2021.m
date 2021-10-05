function [VIF, VIF_adj]=VIF_DNR_2021(X);
% This function calculates the Variance Inflation Factor
%
% USAGE :
% [VIF, VIF_adj]=VIF_DNR_2020(X);
%
% INPUT :
% X = the scaled predictor block (rows, cols)
%
% OUTPUT:
% VIF = vector of Variance Inflation Factors for the variables of X
% VIF_adj = vector of Adjusted Variance Inflation Factors for the variables of X
%
% Calculated doing PLS2_DNR_2020 with (at most) max_LVs
% between the matrix less each variable and the variable removed

% Centre X
Options.CN='N';

%%%%%% ADDED BY DNR
release=version('-release');
if str2num(release(1,1:4))>2014
    parallel='T';
else
    parallel='F';
end
%%%%%% ADDED BY DNR

%%%%%% ADDED BY DNR
if parallel=='T'
    %%%%%% ADDED BY MdF
    % Opens the parallel computing pool if not already opened
    p=gcp('nocreate');
    if isempty(p)==1 % Checks if a pool already exists
        c = parcluster('local');
        c.NumWorkers = feature('numcores'); % No. of workers to be used
        parpool(c, c.NumWorkers);
    end
end

[rows, cols] = size(X);
LVs=1;

%%%%%% ADDED BY DNR
if parallel=='T'
    % SWITCHED TO PARFOR LOOP
    parfor i = 1:cols
        %%%%% CHANGED BY MdF
        [PLS2_Res]=PLS2_DNR_2020(X(:,1:end~=i),X(:,i),LVs,Options); % Just 1 LV
        %%%%% CHANGED BY MdF
        
        y_mean=mean(X(:,i));
        SStot=sum((X(:,i)-y_mean).*(X(:,i)-y_mean));
        SSres=sum((X(:,i)-PLS2_Res.Yhat).*(X(:,i)-PLS2_Res.Yhat));
        R2=1-(SSres./SStot);
        VIF(:,i)=1/(1-R2);
        
        R2_adj=1-((1-R2)*(rows-1)/(rows-LVs-1));
        VIF_adj(:,i)=1/(1-R2_adj);
    end
    
else
    for i = 1:cols
        %%%%% CHANGED BY MdF
        [PLS2_Res]=PLS2_DNR_2020(X(:,1:end~=i),X(:,i),LVs,Options); %
        %%%%% CHANGED BY MdF
        
        y_mean=mean(X(:,i));
        SStot=sum((X(:,i)-y_mean).*(X(:,i)-y_mean));
        SSres=sum((X(:,i)-PLS2_Res.Yhat).*(X(:,i)-PLS2_Res.Yhat));
        R2=1-(SSres./SStot);
        VIF(:,i)=1/(1-R2);
        
        R2_adj=1-((1-R2)*(rows-1)/(rows-LVs-1));
        VIF_adj(:,i)=1/(1-R2_adj);
    end
end
