function [X_CentStd, X_Cent, stdX, mX] = ColCenterStdCal_MZ(X0c);
% For calibration set

[n0c,p]=size(X0c);

% Column center
mX = mean(X0c);

X_Cent = X0c - ones(n0c,1) * mX;

% Column standardize
stdX = std(X_Cent);

X_CentStd= X_Cent ./ (ones(n0c,1) * stdX);


