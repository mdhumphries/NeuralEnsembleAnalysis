function f = bimodalGammaCDF(xdata,a1,b1,a2,b2,p)

% BIMODALGAMMACDF explict form of bimodal Gamma for custom MLE
% F = BIMODALGAMMACDF(X,A1,B1,A2,B2,P1) evaluates the bimodal Gamma
% CDF given data in X and parameters A1,A2 (shape), B1, B2 (scale) and P
% the mixing parameter for mode 1 (P2 = 1 - P1)
%
%
% Mark Humphries 24/4/2013

f1 = gamcdf(xdata,a1,b1);
f2 = gamcdf(xdata,a2,b2);
f = p*f1 + (1-p)*f2;

f(f<=0) = eps;  % handle machine-precision problems when using custom MLE

