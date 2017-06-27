function [bFilt,aFilt,filterFcHz,filterBwHz] = cochlea_linearfilterbank ...
    (Fe, fMinHz, fMaxHz, nFilters, filtOrder, showResults)
% [bFilt,aFilt,filterFcHz,filterBwHz] = COCHLEA_LINEARFILTERBANK 
%           (Fe, fMinHz, fMaxHz, nFilter, filtOrder, showresults)
% Computes the filter coefficients of the cochlea's filterbank. The 
% filterbank is an array of bandpass filter. The center frequencies are
% linearly spaced between fMinHz and fMaxHz.
% The filters used are Butterworth filters. If any unstability is detected,
% the script gives a warning the command line interface.
%
% INPUTS:
%   - Fe            : Sampling frequency (Hz)
%   - fMinHz        : Minimal frequency of analysis for the cochlea (Hz)
%   - fMaxHz        : Maximal frequency of analysis for the cochlea (Hz)
%   - nFilters      : Number of bandpass filters in the filterbank
%   - filtOrder     : Order of these filters (Butterworth filters)
%                     (must be multiple of 2)
%   - showResults   : If 1, show the frequency response of the filterbank
%
% OUTPUTS:
%   - bFilt         : Filter numerator coefficients
%   - aFilt         : Filter denominator coefficients
%   - filterFcHz    : Filter center frequencies (Hz)
%   - filterBwHz    : Bandwidth of each filter (Hz) 
%
% See also: cochlea_erbfilterbank
% 
% Author(s) : Martin Deudon (2016)

if nargin==5
    showResults = 0;
end


%% Filter Design
%- Filter design
filterBwHz      = 2*(fMaxHz-fMinHz)/nFilters;
fBandStart      = linspace(fMinHz,fMaxHz-filterBwHz,nFilters);
fBandEnd        = linspace(fMinHz+filterBwHz,fMaxHz,nFilters);
filterFcHz      = linspace(fMinHz+filterBwHz/2,fMaxHz-filterBwHz/2,nFilters);
fBandStart(fBandStart<0)  = 0;
fBandEnd  (fBandEnd>Fe/2) = Fe/2.01;
aFilt       = zeros(filtOrder+1,nFilters);
bFilt       = zeros(filtOrder+1,nFilters);
for i=1:nFilters
    [bFilt(:,i),aFilt(:,i)] = butter(filtOrder/2,(2/Fe)*[fBandStart(i),fBandEnd(i)],'bandpass');
    %- Check stability
    if ~isstable(bFilt(:,i),aFilt(:,i))
        warning(['Unstable filter detected for index ',num2str(i)]);
    end
end


%% Results
if showResults
    filterInd = 1:nFilters;
    figure;
    plot(filterInd,filterFcHz,'x'); hold on;
    plot(repmat(filterInd,2,1),[fBandStart;fBandEnd],'color','r');
    xlabel('Filter Index'); ylabel('Caracteristic Frequency (Hz)');
    title('Caracteristic frequency of bandpass filters (CCIs) with ERB');

    % Filter impulse responses
    nPoints = 1000;
    H       = zeros(nFilters,nPoints);
    F       = zeros(1,nPoints);
    s.xunits = 'Hz';
    s.yunits = 'db';
    s.plot   = 'mag';
    s.fvflag = 0;
    s.yphase = 'degrees';
    for i=1:nFilters
        [H(i,:),F] = freqz(bFilt(:,i),aFilt(:,i),nPoints,Fe);
    end
    HdB = 20*log10(abs(H));

    figure;
    ax(1) = subplot(211);
    plot(F,HdB);
    axis([fMinHz,Fe/2,-60,5]);
    xlabel('Frequency (Hz)'); ylabel('Gain (dB)'); 
    title(['Filter Bank Impulse Response - Filter order: ',num2str(filtOrder)]);
    ax(2) = subplot(212);
    semilogx(F,HdB);
    axis([fMinHz,Fe/2,-60,5]);
    xlabel('Frequency (Hz)'); ylabel('Gain (dB)'); 
    title(['Filter Bank Impulse Response - Filter order: ',num2str(filtOrder)]);
    linkaxes(ax,'x');
end

end

