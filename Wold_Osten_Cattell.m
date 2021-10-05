function [Wolds_R, Ostens_F, Cattells_RPV] = Wold_Osten_Cattell(PRESS, EigVal, Max_nLVs, nR_Test);
% USAGE :
% [Wolds_R, Ostens_F, Cattells_RPV] = Wold_Osten_Cattell(PRESS, Eig, Max_nLVs, nR_Test);

% INPUT :
% PRESS : Matrix of Prediction Error Sum of Squares
% Eig : Matrix of Eigenvales
% Max_nLVs : Maximum number of LVs (Integer )
% nR_Test : Number of samples (rows) in the Test set

% OUTPUT :
% Wolds_R
% Wold’s R is a vector of the ratios of successive PRESS values.
% The usual cutoff for Wold’s R criterion is when R is greater than unity.
%
% Ostens_F
% A modification of Wold's R

% Cattells_RPV
% Residual Percent Variance (RPV) 
% Assumes that the residual variance should level off, 
% after a suitable number of factors have been extracted
%
% REFERENCES :
%
% *** Wolds_R
% S. Wold
% Cross-validation estimation of the number of components
% in factor and principal component analysis,
% Technometrics 24 (1978) 397– 405.
%
% Model selection for partial least squares regression
% Baibing Li, Julian Morris, Elaine B. Martin
% Chemometrics and Intelligent Laboratory Systems 64 (2002) 79– 89

% *** Ostens_F
% D.W. Osten
% Selection of optimal regression models via crossvalidation,
% J. Chemom. 2 (1988) 39– 48.

% *** Cattells_RPV
% Residual percent variance (RPV; scree test)
% R.B. Cattell
% Multivariate Behav. Res. 1 (1966) 245– 276.

Wolds_R=PRESS(2:Max_nLVs,:)./PRESS(1:Max_nLVs-1,:);
Wolds_R=[Wolds_R;Wolds_R(end,:)];

Ostens_Fup=PRESS(1:Max_nLVs-1,:)-PRESS(2:Max_nLVs,:);
Ostens_Fdn=(PRESS(2:Max_nLVs,:))./(-1*([2:Max_nLVs]'-nR_Test));
Ostens_F=Ostens_Fup./Ostens_Fdn;
Ostens_F=[Ostens_F;Ostens_F(end,:)];

Cattells_RPVup=cumsum(EigVal);
Cattells_RPVdn=sum(EigVal);
Cattells_RPV=Cattells_RPVup./Cattells_RPVdn;
Cattells_RPV=Cattells_RPV(1:Max_nLVs,:);

