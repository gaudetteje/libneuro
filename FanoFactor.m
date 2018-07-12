function F = FanoFactor(S,fs,wLen)
% FANOFACTOR  computes the Fano Factor of a discrete spike train
%
% FANOFACTOR(S,FS,WLENGTH) finds the spike count <n> for a spike train, S,
%  based on non-overlapping windows of duration wLen.  The Fano Factor is
%  then found using the formula F = std(n)/mean(n).

% force S to be a column vector
S(:)=S;

% reshape spike train into sequential matrix of specified duration
nSamples = ceil(wLen*fs);
nWindows = floor(numel(S)/nSamples);
N = reshape(S(1:nSamples*nWindows),nSamples,nWindows);  % truncate remaining samples

N = sum(N,1);   % sum along time axis for each windowed sequence

F = var(N) / mean(N);
