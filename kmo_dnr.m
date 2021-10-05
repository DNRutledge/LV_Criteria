function [kmo] = kmo_dnr(X)
%KMO Kaiser-Meyer-Olkin Measure of Sampling Adequacy.
% Sampling adequacy predicts if data are likely to factor well,
% based on correlation and partial correlation.
% It has been suggested that inv(R) should be a near-diagonal matrix in order
% to successfully fit a factor analysis model.
% To assess how close inv(R)
% is to a diagonal matrix, Kaiser (1970) proposed a measure of sampling
% adequacy, now called KMO (Kaiser-Meyer-Olkin) index. The common part, called
% the image of a variable, is defined as that part which is predictable by
% regressing each variable on all other variables.
%
% A KMO index <= 0.5 indicates the correlation matrix
% is not suitable for factor analysis.
%
% USAGE :
% [kmo] = kmo_dnr(X);
%
% INPUT :
% X - Input matrix can be a data matrix (size n-data x p-variables)
%
% OUTPUT :
% kmo :  - Kaiser-Meyer-Olkin Index.
%
% REFERENCE :
% Rencher, A. C. (2002)
% Methods of Multivariate Analysis. 2nd. ed.
% New-Jersey:John Wiley & Sons. Chapter 13 (pp. 408-450).
%

% Changed by DNR 31/12/2017
% error(nargchk(1,1,nargin));
error(nargchk(1,3,nargin));

% Changed by DNR 31/12/2017
% msg = nargoutchk(1, 2, nargout);
msg = nargoutchk(1, 3, nargout);

% Added by DNR 31/12/2017
% A & B are NaNs if there is a variable which is zero for all individuals
% So replace the zeros by very small random values
Vars_zero=find(sum(X)==0);
if Vars_zero>0
    X_temp=Zero2Randn(X, 100000);
    X(:,Vars_zero)=X_temp(:,Vars_zero);
end

X = corrcoef(X);

% Changed by DNR 31/12/2017
iX = pinv(X);
% iX = pinv_DNR(X);
% % % iX = inv(X);

S2 = diag(diag((iX.^-1)));
AIS = S2*iX*S2; %anti-image covariance matrix

% But what if AIS has negative values ?
Dai = diag(diag(sqrt(AIS)));

%anti-image correlation matrix
AIR = (Dai)\AIS/(Dai);

a = sum((AIR - diag(diag(AIR))).^2);
AA = sum(a);
b = sum((X - eye(size(X))).^2);
BB = sum(b);

N = BB;
D = AA+BB;
kmo = N/D;


return;