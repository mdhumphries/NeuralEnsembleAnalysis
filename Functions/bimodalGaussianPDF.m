function f = bimodalGaussianPDF(xdata,mu1,std1,mu2,std2,p)

% BIMODALGAUSSIANPDF explict form of bimodal Gaussian for custom MLE
% F = BIMODALGAUSSIANPDF(X,M1,S1,M2,S2,P1) evaluates the bimodal Gaussian
% PDF given data in X and parameters M1,M2 (means), S1, S2 (std dev) and P
% the mixing parameter for mode 1 (P2 = 1 - P1)
%
%
% Mark Humphries 23/4/2013

% f1 = 1./(std1.*sqrt(2*pi)) .* exp(-(xdata - mu1).^2/(2*std1.^2));
% f2 = 1./(std2.*sqrt(2*pi)) .* exp(-(xdata - mu2).^2/(2*std2.^2));

f1 = normpdf(xdata,mu1,std1);
f2 = normpdf(xdata,mu2,std2);
f = p*f1 + (1-p)*f2;

f(f<=0) = eps;

if any(isnan(f)) | any(isinf(f))
     keyboard
end