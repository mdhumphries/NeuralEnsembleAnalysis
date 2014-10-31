function bic = BICL(L,N,K)

%BICL Bayesian Information Criterion for MLE model comparison
%   BICL(L,N,K) where L is the negative log-likelihood value for the model.
%   N is the number of data points and K is the number of coefficients. 
%   Computes and returns the BIC score for the model.
%   
%   The model with the lowest BIC value is the best fit. 
%
%   References: 
%   Schwarz, G. (1978) Estimating the dimension of a model Annals of Statistics, 6, 461-464
%
%   Wagenmakers, E.-J. & Farrell, S. (2004) AIC model selection using
%   Akaike weights. Psychon Bull Rev, 2004, 11, 192-196
%
%   NOTES: 
%   (1)this computes BIC from likelihood (L); use BICSS for estimate from
%   residuals (sum square error).
%
%   (2) The penalty for extra parameters is stiffer in BIC than in AIC
%
%   (3) "Schwarz" weights for BIC can be computed in exactly the same way
%   as the Akaike's weights for AIC: simply supply the BIC scores returned
%   by this function to the AKWGTS function (note: could also use BICWGTS function
%   but the AKWGT function scales by score magnitude, so is less likely to return
%   NaNs).
%
%   Mark Humphries 12/05/2011


% raw BIC
% bic = -2 * log(L) + K*log(N);    % where L is likelihood

bic = 2.* L + K.*log(N); % where L is negative log-likelihood [most likely scenario]
