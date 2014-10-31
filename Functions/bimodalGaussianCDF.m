function f = bimodalGaussianCDF(xdata,mu1,std1,mu2,std2,p)

% BIMODALGAUSSIANCDF explict form of bimodal Gaussian CDF for custom MLE
% F = BIMODALGAUSSIANCDF(X,M1,S1,M2,S2,P1) evaluates the bimodal Gaussian
% CDF given data in X and parameters M1,M2 (means), S1, S2 (std dev) and P
% the mixing parameter for mode 1 (P2 = 1 - P1)
%
%
% Mark Humphries 24/4/2013

f1 = normcdf(xdata,mu1,std1);
f2 = normcdf(xdata,mu2,std2);
f = p*f1 + (1-p)*f2;

f(f<=0) = eps; % handle machine-precision problems when using custom MLE

if any(isnan(f)) | any(isinf(f))
     keyboard
end