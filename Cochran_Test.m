function [Cochran] = Cochran_Test(X)
%
% USAGE : 
% [Cochran] = Cochran_Test(X)
%
% INPUT : 
% X - data matrix to be tested for sphericity
%
% OUTPUT :
% Cochran : Matrix with the values for the Cochran hypothesis test
% of ratio of max(Var)/sum(Var) of the variables
%
% REFERENCE :
% W.G. Cochran
% The distribution of the largest of a set of estimated variances as a fraction of their total
% Annals of Human Genetics (London) 11(1), 47–52 (January 1941)

VarX=var(X);
maxVar=max(VarX);
sumVar=sum(VarX);

Cochran=maxVar/sumVar;



