function aic = AICL(L,N,K)

%AICL Akaike's Information Criterion for MLE model comparison
%   AICL(L,N,K) where L is the negative log-Likelihood value for the model.
%   N is the number of data points and K is the number of coefficients. 
%   Computes and returns the corrected AIC score for the model.
%   
%   The model with the lowest AIC value is the best fit.
%
%   References: (1) Burnham, K. P. and Anderson, D. R. (2002) Model Selection and 
%   Multimodel Inference: A Practical-Theoretic Approach, 2nd ed.
%   Springer-Verlag  [TO GET]
%
%   (2) Akaike, H. (1974) "A new look at the statistical model identification." IEEE Transactions on Automatic Control, AC-19, 716-723
%
%   (3) Hurvich, C. M., and Tsai, C-L. (1989). Regression and time series model selection in small samples. Biometrika, 76, 297-307.
%   [the AIC correction]   
%
%   NOTE: this computes AIC from negative log-likelihood (L); least-squares curve-fits can also be compared when 
%   using SSE (sum-square error) using AICSS instead!
%
%   Mark Humphries 12/05/2011

% K  = K + 1; % additional degree-of-freedom is model!

% raw AIC
% aic =  -2 * log(L) + 2 * K;    % where L is likelihood
aic = 2.*L+2.*K;                   % where L is already negative log-likelihood (the likely case)
% apply correction in case N close to K
aic = aic + (2.*K.*(K+1)) ./ (N - K - 1);
