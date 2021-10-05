function [PLS2_Res]=PLS2_DNR_2020(X,Y,lvs, Options);
% This function implements PLS2
%   USAGE:
%        [PLS2_Results]=PLS2_DNR_2020(X,Y,lvs);
%   INPUTS:
%        X = the scaled predictor block (n,px)
%        Y = the scaled predictand block (n,py)
%        lvs = the number of latent variables to be calculated (1,1)
%        Options.CN='N' or 'C' (default) Data Centred or not 
%
%   OUTPUTS:
% PLS2_Results : A structure containing :
%        B = the matrix of regression vectors or matrices (px,lvs)
%        B0 = regression intercepts (1,py)
%        T = X scores (n,lvs)
%        P = X loadings (px,lvs)
%        W = X weights (px,lvs)
%        U = Y scores (n,lvs)
%        Q = Y loadings (py,lvs)
%        Y_hat = Predicted Y values (n, py)

if exist('Options','var')
    if isfield(Options,'CN')
        CN=Options.CN;
    else
        CN='C';
    end
end    

X_old=X;
[rows,px] = size(X);
[rows,py] = size(Y);
u_old = Y(:,1);

W=zeros(px,lvs);
P=zeros(px,lvs);
Q=zeros(py,lvs);
T=zeros(rows,lvs);
U=zeros(rows,lvs);

%%%%% COMMENTED BY DNR
% Variable size because maybe PLS1 or PLS2
% B=zeros(px,py);

meanX=mean(X); %DB, 05.04.16
meanY=mean(Y); %DB, 05.04.16

%  DNR 01/06/17
if CN=='C'
    X = X - ones(rows,1)*meanX;
end

% Necessary for PLS_ICA where nLVs decreases to 1
B_mat=zeros(lvs,px,py);
B0_mat=zeros(lvs,py);

for i = 1:lvs
    limite=0;
    while ( 1 )
        w = u_old' * X;
        w = w /norm(w);
        t = X * w';
        t = t / (w*w');
        q = t' * Y;
        q = q / (t'*t);
        u = (Y * q')/(q*q');
        if ( norm(u - u_old) < 1e-6)
            break;
        end
        limite=limite+1;
        if ( limite > 1000)
            break;
        end
        u_old = u;
    end
    p = t'*X;
    p = p / (t'*t);
    
    W(:,i)=w';
    P(:,i)=p';
    Q(:,i)=q';
    T(:,i)=t;
    U(:,i)=u;
    
    X = X - t*p;
    Y = Y - t*q;

end

B = W/(P'*W)*Q';
B0 = meanY - meanX*B; %DB, 05.04.16

PLS2_Res.B=B;
PLS2_Res.B0=B0;
PLS2_Res.Scores=T;
PLS2_Res.P=P;
PLS2_Res.W=W;
PLS2_Res.U=U;
PLS2_Res.Q=Q;

PLS2_Res.Yhat=X_old*B+ones(size(X_old,1),1)*B0;

