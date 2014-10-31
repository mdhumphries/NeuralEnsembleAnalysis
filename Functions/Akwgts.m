function w = Akwgts(AICs)

% AKWGTS computes Akaike's weights from a set of AIC values
%   W = AKWGTS(A) computes Akaike's weights W, from an array of AIC values
%   given in A. These give the relative likelihood of the model
%
%   Reference:
%   Wagenmakers, E.-J. & Farrell, S. (2004) AIC model selection using Akaike
%   weights. Psychon Bull Rev, 11, 192-196
%
%   Burnham, K. P., and D. R. Anderson. (2002). Model Selection and Multimodel Inference: 
%   a practical information-theoretic approach, 2nd edition. Springer-Verlag, New York.
%
% Mark Humphries 19/4/2012

w = zeros(numel(AICs),1);

% deal with NaNs
ixnans = find(isnan(AICs));
ixNnans = find(~isnan(AICs));
AICs = AICs(ixNnans);
n = numel(AICs);
bestAIC = min(AICs);

deltaAICs = AICs - bestAIC;

% compute weights
w(ixnans) = 0;
w(ixNnans) = exp(-0.5 * deltaAICs) ./ sum(exp(-0.5 * deltaAICs));