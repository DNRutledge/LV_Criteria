function [RE,IE,XE,IND] = Malinowski(Ri, i);
% USAGE :
% [RE,IE,XE,IND] = Malinowski(Ri, i);

% INPUT :
% Ri : Matrix of residuals
% after removing i components

% OUTPUT :
% RE,IE,XE,IND

% REFERENCES :
% Malinoiwsi's criteria
%     E. R. Malinowski
% Theory of Error in Factor Analysis
% Anal. Chem., 49(4), 606-612  (1977).
%
%     E. R. Malinowski
% Determination of the Number of Factors and the Experimental Error in a Data Matrix
% Anal. Chem., 49(4), 606-612-617 (1977).

[nR,nC]=size(Ri);

if nR>nC+1
    nX=nR;
    nY=nC;
else
    nX=nC;
    nY=nR;
end
    

E2=Ri.*Ri;
RSD2=sum(E2(:))/(nX*(nY-i));
RE=sqrt(RSD2);
IE=sqrt(RSD2*i/nY);
XE=sqrt(RSD2*(nY-i)/nY);
IND=sqrt(RSD2/(nY-i));


