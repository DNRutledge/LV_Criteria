function [Hartley] = Hartley_Test(X);
% USAGE : 
% [Hartley] = Hartley_Test(X);
%
% INPUT : 
% X - data matrix to be tested for sphericity
%
% OUTPUT :
% Hartley : Matrix with the values for the Hartley hypothesis test
% of ratio of max(Var)/min(Var) of the variables
%
% REFERENCE :
% Hartley, H.O. (1950).
% The maximum F-ratio as a short cut test for homogeneity of variance,
% Biometrika, 37, 308-312.
%
VarX=var(X);
maxVar=max(VarX);
minVar=min(VarX);

Hartley=maxVar/minVar;



