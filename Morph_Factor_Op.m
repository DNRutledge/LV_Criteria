function [MFactor] = Morph_Factor_Op(dir,X);
% [MFactor] = Morph_Factor_Op(dir,X);
%
% Calculates the morphological Factor for a matrix X.
%
% 23.04.2020
%
% INPUT :
% dir   'row' for calculation by rows and any other string or value for
%       calculation by columns
% X     signal (e.g., a loading vector matrix, or a matrix of vectors of
%       regression coefficients)
%
% OUTPUT :
% MFactor : Matrix with the values for the Morphological Factor
% for each row or column

% REFERENCES :
% Chemical rank estimation by multiresolution analysis for two-way data
% in the presence of background
% Hailin Shen, Jihong Wang, Yizeng Liang, Karin Pettersson, Mats Josefson, Johan Gottfries, Frank Lee
% Chemometrics and Intelligent Laboratory Systems 37(2) 1997 261-269
%
% The morphological score and its application to chemical rank determination
% Hailin Shen, Laila Stordrange, Rolf Manne, Olav M. Kvalheim, Yizeng Liang
% Chemometrics and Intelligent Laboratory Systems 51 2000 37–47

if strcmp(dir,'row')
    % By Rows
else
    % By Columns
    X=X';
end

[n, p] = size(X);

MeanA = mean(X');
A_MeanA =X-MeanA'*ones(1,p);

for i=1:n
    NormA_meanA=norm(A_MeanA(i,:));
    
    DiffA_MeanA=diff(A_MeanA(i,:));
    pos = DiffA_MeanA>0;
    changes = xor(pos(1:end-1),pos(2:end));
    % zero crossing points 
    ZCP = sum(changes); 
    
    NormDiffA_MeanA=norm(DiffA_MeanA);
    
    MFactor(i,:)=NormA_meanA/NormDiffA_MeanA/ZCP;
end

