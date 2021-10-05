function [Bartlett_X2] = Bartlett_DNR(X);
%
% USAGE : 
% [Bartlett_X2, Proba] = Bartlett_DNR(X);
%
% INPUT : 
% X - data matrix to be tested for sphericity
%
% OUTPUT :
% Bartlett_X2 : Bartlett's X²
%
% REFERENCE :
% Bartlett, M. S. (1951).
% The Effect of Standardization on a X² Approximation in Factor Analysis.
% Biometrika 38(3/4), 337–344.
%
%  Bartlett’s test for sphericity
% (testing that the correlation matrix has an identity matrix)
% Bartlett’s test for Sphericity compares your correlation matrix
%(a matrix of Pearson correlations)
% to the identity matrix.
% In other words, it checks if there is a redundancy between variables
% that can be summarized with some factors
%
% X² = -[(n-1)-(2k + 5)/6] * log(|R|) 
% where :
% n is the number of observations,
% k the number of variables,
% and R the correlation matrix of the data supplied in X.
% |R| is the determinant of R.
% 
% Bartlett's X² is asymptotically X²-distributed
% with df = k*(k-1)/2 under the null hypothesis.
% Note that, because the bias-corrected correlation matrix is used,
% (n-1) is employed instead of n, as in the paper.
%
%%%%%% DNR
% I used the original formulae for df & X²
%%%%%% DNR
[nR,nC]=size(X);

% R=corrcoef(X);
% det_R=det(R);

% Bartlett_X2=-((nR-1)-(2*nC+5)/6)*log(det_R);
% Negative determinants ! with old R !
% Bartlett_X2=-((nR-1)-(2*nC+5)/6)*log(abs(det_R));


[X_cs, X_Cent, stdX, mX] = ColCenterStdCal_MZ(X);
R=X_cs'*X_cs;

det_R=det(R);

% % % % % Bartlett_X2=-(nR-(2*nC+5)/6)*log(det_R);
% Negative determinants ! with old R !
Bartlett_X2=-((nR-1)-(2*nC+5)/6)*log(abs(det_R));


% Correct for df ?
% df=((nC+2)*(nC-1)/2); %degrees of freeedom
df=(nC*(nC-1)/2); %degrees of freeedom

Bartlett_X2=Bartlett_X2/df;

% Proba = 1-chi2cdf(Bartlett_X2,df);  %Probability that null Ho: is true.


