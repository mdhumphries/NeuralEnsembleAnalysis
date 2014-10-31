function [y,bins,f1,f2] = LIF_xcorr(times1,times2,bin_size,time_window,varargin)

% LIF_XCORR cross-correlogram/covariogram histogram
%
%   LIF_XCORR(A,B,BINSIZE,T) where A and B are time-stamp (in seconds) arrays of spiking events,
%   BINSIZE is the size of the sample bin in seconds, T is a 2-element array specifying the start and 
%   end of the spike trains in seconds. Computes the cross-correlogram between two spike trains, 
%   with A as the reference train. This is most useful when the spike trains can be compared over their entirety. 
%   For computing correlations between sets of spike trains recorded over many trials, use either
%   JPSTH or COVARIOGRAM (the latter calls this function)
%
%   LIF_XCORR(...,MAXLAG) sets the size of the maximum 
%   interval returned to be +/- MAXLAG seconds; 
%
%   [Y,BINS,F1,F2] = LIF_XCORR(...) returns the binned interval counts Y and the bin centres B. 
%   Plot using BAR(BINS,Y,1). Can also optionally specify F1 and F2 to return the mean firing rates of trains A and B.
%   
%   NOTE#1: specify a maximum interval MAXLAG for best results. It ensures that the histogram is
%   approximately flat if the signals are not correlated. If the whole
%   signal is used for both reference and comparison, then edge effects
%   will result as the extreme intervals are necessarily low in frequency.
%
%   NOTE#2: try initial values of BINSIZE = 0.001 and MAXLAG = 0.1 and
%   adjust as necessary - BINSIZE should generally be the quantising step of the underlying spike-trains,
%   or a sufficiently small time-window to guarantee a single spike per bin in each train; so for simulated trains,
%   this could be the smallest absolute refractory period if known
%
%   NOTE#3: to get autocorrelogram, just enter the same spike-train for A and B
%
%   REFERENCE: Dayan, P & Abbott, L. F. (2001) Theoretical Neuroscience. Cambridge, MA: MIT Press
%
%   Mark Humphries 6/10/07

if time_window(2) < time_window(1)
    error('End of spike train time window must be after start')
end

time_seconds = time_window(2) - time_window(1);

% turn arrays right way round for further processing
[r1 c1] = size(times1);
[r2 c2] = size(times2);
if r1 > c1
    times1 = times1';
end
if r2 > c2
    times2 = times2';
end

if nargin >= 5 & isnumeric(varargin{1}) & ~isempty(varargin{1})
     max_lag = varargin{1};
else
    max_lag = time_seconds;    % do not use max lag
end

% bins = -time_seconds:bin_size:time_seconds;
bins = -max_lag:bin_size:max_lag; 

num_spikes1 = length(times1);
num_spikes2 = length(times2);

f1 = num_spikes1/time_seconds;
f2 = num_spikes2/time_seconds;

if num_spikes1 < 2 | num_spikes2 < 2
    % then not sufficient to process
    y = zeros(1,length(bins));
    return
end

%%%% find difference between spike times....
mat1 = repmat(times1,length(times2),1);
mat2 = repmat(times2',1,length(times1));

% keyboard

diffT = mat2 - mat1;             % matrix of all time differences between spikes
% y = hist(diffT(:),bins);        % bin time-differences according to bin array
y = hist(diffT(abs(diffT)<=max_lag),bins);        % bin time-differences according to bin array

%keyboard

% y = y - (num_spikes1 * num_spikes2)/numel(bins);
% % ycov = y - (f1*f2);
% % ynorm = y ./ ((time_seconds - abs(bins)) .* sqrt(f1*f2));

% % reduce to requested range: ditch first and last bins for "garbage"
% if max_lag ~= time_seconds
%      bins = bins(2:end-1);  
%      y = y(2:end-1);
% end









