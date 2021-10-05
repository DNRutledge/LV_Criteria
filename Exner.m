function [Exners_PHI] = Exner(R1, Ri, i);
% USAGE :
% [Exners_PHI] = Exner(R1, Ri, i);

% INPUT :
% R1 : Original undeflated Matrix
% Ri : Matrix of residuals
% after removing i components

% OUTPUT :
% Exners_PHI

% REFERENCES :
% Exner, O. (1966).
% Additive physical properties. I.
% General relationships and problems of statistical nature.
% Collection of Czech, Chem. Commun. 31, 3222-3251. doi: 10.1135/cccc19663222

[nR,nC]=size(Ri);

if nR>nC+1
    nX=nR;
    nY=nC;
else
    nX=nC;
    nY=nR;
end
    

% Exner function (PHI).
GrandMean=mean(R1(:));
Exnerup=nR*nC*sum(sum(Ri.*Ri));
Exnerdn=(nR*nC-i)*sum(sum((R1-GrandMean).*(R1-GrandMean)));
Exners_PHI=Exnerup./Exnerdn;

    

