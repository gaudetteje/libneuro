function [H,M] = autocorrspike(S,fs,wLen)
% AUTOCORRSPIKE  generates the spike-train autocovariance histogram
%
% [H,M]=AUTOCORRSPIKE(SPIKES,FS,WINLENGTH,BINWIDTH) returns the autocorrelation function
% of the neural response function with its average over time/trials removed
%







% calculate window length in samples
nPoints = ceil(wLen*fs);

% compute autocorrelation function for specified window length
H = xcorr(S,nPoints);

% create histogram bins
M = [-wLen : 1/fs : wLen]';

% normalize bins and subtract the mean of a uniformly distributed spike train
n = sum(S);
T = length(S)/fs;
H = H./T - ((n^2 / fs) / T^2);
H(nPoints+1) = 0;       % remove peak at 0

% % plot histogram if no output arguments exist
% if ~nargout
%     bar(H,M)
% end

%% we can use the reshape function based on how many points are within bins
%% of a specified width...  Then plot the bar graph of the resulting
%% vector.