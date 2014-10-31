function [coeffs,AICs,BICs,pAICs,pBICs,L] = fitMLEdistribution(xdata,fits,varargin)

% FITMLEDISTRIBUTION find MLE fits of all common distributions, and compare
%   [C,A,B,pA,pB,L] = FITMLEDISTRIBUTION(X,FITS) takes the univariate data in vector X, 
%   and finds maximum likelihood estimators (using the 'MLE' functions from
%   Statistics Toolbox) for all the distributions specified in array FITS. It then
%   computes AIC and BIC scores for all fits, and the corresponding estimates of p(model).
%   The set of distributions is defined by the array FITS, whose numeric entries 
%   request the following:
%       1 = normal
%       2 = lognormal
%       3 = gamma
%       4 = exponential
%       5 = Weibull
%       6 = power-law   (requires extra parameter [X,FITS,M], where M is
%       the minimum value of X for which the power-law is thought to hold)
%       7 = uniform
%       8 = bimodal normal [returns as (mu1, std1, mu2, std2, p)]
%       9 = bimodal gamma [returns as (a1,b1,a2,b2,p)]   
%
%   Returns: C a cell array of the coefficients of all the distribution;
%   A,B the arrays of AIC, BIC scores for  every tested distribution; pA, pB - the esimates of p(model) 
%   based respectively on the AIC and BIC score; L the negative log-likelihood value for each 
%   distribution; 
%
%   Notes:
%   (1)  Other distributions can be added (for example, truncated power-law distributions can be
%   added from Khanin & Wit (2006) with a little work)
%
%   (2) The power law MLE solution is from Newman, M. E. J. (2005) Power laws,
%   Pareto distributions and Zipf’s law. Contemporary Physics, 46, 323-351.
%
%   (3) Bimodal Gaussian and Gamma use the custom options from the MLE
%   function, and specify "FMINCON" as the search function - this requires
%   the Optimization Toolbox.
%
%   (4) Approximate p(model) scores for both AIC and BIC are computed using
%   the AKWGTS function, which is based on the size of the score
%   differences rather than on the absolute values of the score
%
%   Mark Humphries 22/4/2013

if any(fits > 9)
    error('Non-existent distribution chosen')
end

if find(fits == 6) & nargin < 3
    error('Please specify parameter M for power-law distribution fit')
end

nfits = length(fits);
coeffs = cell(nfits,1);
L = zeros(nfits,1);

% options for custom MLE in case we need them
options = statset('mlecustom');
% options.MaxIter = 1000; options.MaxFunEvals = 1000;

for loop = 1:nfits
    switch fits(loop)
        case 1
            % normal
            [m,s] = normfit(xdata);
            coeffs{loop} = [m,s];
            Y =  pdf('normal',xdata,m,s);
            L(loop) = -sum(log(Y));

        case 2
            % log-normal
            coeffs{loop} = lognfit(xdata);
            Y =  pdf('lognormal',xdata,coeffs{loop}(1),coeffs{loop}(2));
            L(loop) = -sum(log(Y));
        case 3
            % gamma
            coeffs{loop} = gamfit(xdata);
            Y =  pdf('gamma',xdata,coeffs{loop}(1),coeffs{loop}(2));
            L(loop) = -sum(log(Y));

        case 4
            % exponential
            coeffs{loop} = expfit(xdata);
            Y =  pdf('exponential',xdata,coeffs{loop}(1));
            L(loop) = -sum(log(Y));
 
        case 5
            % Weibull
            coeffs{loop} = wblfit(xdata);
            Y =  pdf('weibull',xdata,coeffs{loop}(1),coeffs{loop}(2));
            L(loop) = -sum(log(Y));

        case 6
            % power-law
            coeffs{loop} = pwrfit(xdata,varargin{1});
            L(loop) = pwrlike(coeffs{loop},xdata,varargin{1});
        case 7
            % uniform (2 parameter)
            coeffs{loop} = mle(xdata,'distribution','uniform');
            Y = pdf('uniform',xdata,coeffs{loop}(1),coeffs{loop}(2));
            L(loop) = -sum(log(Y));

        case 8
            try
            % bimodal normal
                coeffs{loop} = mle(xdata,'pdf',@bimodalGaussianPDF,'cdf',@bimodalGaussianCDF,...
                    'start',[mean(xdata)/2,std(xdata)/2,mean(xdata)*2,std(xdata),0.5],...
                    'lowerbound',[0 1e-4 0 1e-4 0],'upperbound',[inf inf inf inf 1],...
                    'options',options,'optimfun','fmincon');
                Y = bimodalGaussianPDF(xdata,coeffs{loop}(1),coeffs{loop}(2),coeffs{loop}(3),coeffs{loop}(4),coeffs{loop}(5));
                L(loop) = -sum(log(Y));
            catch
                % can hard-fail during algorithm
                coeffs{loop} = nan;
                L(loop) = nan;
                keyboard
            end
        case 9
            % bimodal Gamma distribution
            try
                coeffs{loop} = mle(xdata,'pdf',@bimodalGammaPDF,'cdf',@bimodalGammaCDF,...
                    'start',[mean(xdata),std(xdata),max(xdata),1,0.5],...
                    'lowerbound',[1e-4 1e-4 1e-4 1e-4 0],'upperbound',[inf inf inf inf 1],...
                    'options',options,'optimfun','fmincon');
                Y = bimodalGammaPDF(xdata,coeffs{loop}(1),coeffs{loop}(2),coeffs{loop}(3),coeffs{loop}(4),coeffs{loop}(5));
                L(loop) = -sum(log(Y));
            catch
                keyboard
            end
    end

end


%%% do AIC 
AICs = zeros(1,nfits);
for j = 1:nfits
    AICs(j) = AICL(L(j),numel(xdata),numel(coeffs{j}));
end
pAICs = Akwgts(AICs);

%%% do BIC
BICs = zeros(1,nfits);
for j = 1:nfits
    BICs(j) = BICL(L(j),numel(xdata),numel(coeffs{j}));
end
pBICs = Akwgts(BICs);




