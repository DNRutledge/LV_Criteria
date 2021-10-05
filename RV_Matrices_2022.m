function [RV, COI, RIP] = RV_Matrices_2022(X1,X2);
% USAGE :
% [RV, COI, RIP] = RV_Matrices_2022(X1,X2);
% This function calculates 4 matrix correlationvalues between X1 & X2
%
% INPUT:
%        X1 = Matrix of size (n, p)
%        X2 = Matrix of size (n, q)
%
% OUTPUT:
% The RV, COI and RIP coefficients reflect the similarity of the spatial configurations
% of the individuals associated with X1 and X2
%
% ***** RV
% % It ranges between 0 and 1.
% % The lower bound is reached if all the variables in X are uncorrelated with those in X2
% % RV = 1 if and only if X1 and X2 can be matched by a rotation
% % and a multiplication by a scaling factor.

% REFERENCES :
% Stéphane Dray
% "On the number of principal components: A test of dimensionality
% based on measurements of similarity between matrices"
% Computational Statistics & Data Analysis 52 (2008) 2228 – 2237
%
% Robert, P.; Escoufier, Y. (1976). "A Unifying Tool for Linear
% Multivariate Statistical Methods: The RV-Coefficient".
% Applied Statistics. 25 (3): 257–265

% ***** COI
% According to S. Dray
% The numerator of the RV coefficient corresponds to the co-inertia criterion
% which is a measurement of the link between the two tables

% REFERENCES :
% Stéphane Dray, Daniel Chessel & Jean Thioulouse (2003)
% Procrustean co-inertia analysis for the linking of multivariate datasets,
% Écoscience, 10:1, 110-119,
% DOI: 10.1080/11956860.2003.11682757

% ***** RIP
% Inner Product Correlation
% It ranges between 0 and 1

% REFERENCES :
% Ramsay, J.O., J.M.F. Ten Berge and G.P.H. Styan,
% Matrix correlation, Psychometrika, 49 (1984) 403-423.
%
% and studied by Kiers et al. 
% Kiers, H.A.L, Cléroux, R. & Ten Berge, M.F. (1994)
% Generalized analysis based on optimizing matrix correlations and a relation with IDIOSCAL.
% Computational Statistics and Data Analysis : 18, 331-340.

% ***** RLS
% % It should range between 0 and 1 !

% REFERENCES :
% The RLS index is attributed to Lingoes & Schönemann
% Lingoes, J.C. & Schönemann, P.H. (1974)
% Alternative measures of fit for the Schönemann-Carrol matrix fitting algorithm.
% Psychometrika : 39, 423-427.
%
% by Lazraq & Coll.
% Lazraq, A., Cléroux, R. & Kiers, H.A.L. (1992)
% Mesures de liaison vectorielle et généralisation de l'analyse canonique.
% Revue de Statistique Appliquée : 39, 23-35
%
% and studied by Kiers et al. 
% Kiers, H.A.L, Cléroux, R. & Ten Berge, M.F. (1994)
% Generalized analysis based on optimizing matrix correlations and a relation with IDIOSCAL.
% Computational Statistics and Data Analysis : 18, 331-340.


%%
[n,p]=size(X1);
[n,q]=size(X2);
% p==q !

% [X1, Norm, mX]=Normalise_DB(X1); % MS
% [X2, Norm, mX]=Normalise_DB(X2); % MS
% 
[X1,mX] = ColMeanCenterCal_MZ(X1);
[X2,mX] = ColMeanCenterCal_MZ(X2);

% According to S. Dray
% COI=trace((X1'*X2)*(X2'*X1));

% According to S. Dray
% RV=trace((X1*X2')*(X2*X2'))/(sqrt(trace((X1'*X1)*(X1'*X1))*trace((X2'*X2)*(X2'*X2))));
% RV=COI/(sqrt(trace((X1'*X1)*(X1'*X1))*trace((X2'*X2)*(X2'*X2))));

% According to D. Chessel
COI=trace((X1*X1')*(X2*X2'));
RV=COI/(sqrt(trace((X1*X1')*(X1*X1'))*trace((X2*X2')*(X2*X2'))));

% According to S. Dray
% And according to Kiers
% RLS=trace(sqrt((X1'*X2)*(X2'*X1)))/sqrt(trace(X1'*X1)*trace(X2'*X2));

% According to D. Chessel
% RLS=trace(sqrt((X1*X1')*(X2*X2')))/sqrt(trace(X1*X1')*trace(X2*X2'));

% Ramsay, J.O., J.M.F. Ten Berge and G.P.H. Styan,
% Matrix correlation, Psychometrika, 49 (1984) 403-423.
RIP=trace(X1'*X2)/sqrt(trace(X1'*X1)*trace(X2'*X2));




