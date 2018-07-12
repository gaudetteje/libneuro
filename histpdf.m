function varargout = histpdf(data,varargin)
% HISTPDF  plots the histogram normalized to approximate a PDF (area integrates to 1)
%
% HISTPDF(DATA) plots the normalized histogram (area = 1), using 50 bins
%     spread linearly across the min/max data points.
% HISTPDF(DATA,BINS) uses the number of bins specified in BINS
% HISTPDF(DATA,BINS,RANGE) accepts a 2 element vector, RANGE, specifying
%     the desired min/max histogram limits.
% [H,E] = HISTPDF(DATA) returns the normalized histogram values, H, and
%     edge locations, E.
%
% Note:  Bin width can be calculated using BINRANGE/NBINS.

switch nargin
    case 1
        nBins = 50;
        bRange = [min(data) max(data)];
    case 2
        nBins = varargin{1};
        bRange = [min(data) max(data)];
    case 3
        nBins = varargin{1};
        bRange = varargin{2};
    otherwise
        error('HISTPDF requires 1, 2, or 3 input parameters.  Try "help histpdf" for usage.')
end

% check for badness
if min(data) < bRange(1)
    nMax = length(data(data < bRange(1)));
    warning(sprintf('%d data points are less than the specified histogram range and have been excluded.',nMax))
end

if max(data) > bRange(2)
    nMax = length(data(data > bRange(2)));
    warning(sprintf('%d data points are greater than the specified histogram range and have been excluded.',nMax))
end

% create vector of bin edges and run through 'histc'
bEdges = linspace(bRange(1),bRange(2),nBins+1);
h = histc(data,bEdges);

% always normalize histogram to have area = 1
h_norm = diff(bRange)*numel(data)/nBins;

% plots to new or current figure if no return values
if ~nargout
        bar(bEdges, h./h_norm, 'histc');
        set(gca,'XLim',bRange);
        grid on;
else
        varargout(1) = {h./h_norm};
        varargout(2) = {bEdges};
end
