function [DW] = DurbWatMatrix_Op(dir,X);
%
% USAGE :
% [DW] = DurbWatMatrix_Op(dir,X);
% 
% 5.7.2011
% 
% Calculates the Durbin-Watson criterion for a matrix X.
%
% INPUT
% dir   'row' for calculation by rows and any other string or value for 
%       calculation by columns
% X     signal (e.g., a loading vector matrix, or a matrix of vectors of 
%       regression coefficients)
%
% OUTPUT
% DW row vector whith the values for the Durbin-Watson criterion for each
% row or column
%           dw = 0 : there is a strong correlation between
%           successive points
%           dw = 2 : there is a weak correlation (random distribution)
%           between successive points.
%
% REFERENCES :
% Durbin J, Watson GS. (1971)
% Testing for serial correlation in least squares regression. III.
% Biometrika, 37, 1-19., doi: 10.2307/2334313
%
% Rutledge, D. N., Barros A. S. (2002)
% The Durbin-Watson statistic as a morphological estimator of information content.
% Analytica Chimica Acta 446, 279-294. doi: 10.1016/S0003-2670(01)01555-0
%
if strcmp(dir,'row')
   B=((X(:,2:end)-X(:,1:end-1)).^2);
   C=(X.^2);
   DW=(sum(B,2)./sum(C,2))';
   
else
   B=((X(2:end,:)-X(1:end-1,:)).^2);
   C=(X.^2);
   DW=sum(B,1)./sum(C,1);
end
