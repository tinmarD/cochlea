function [bFilt,aFilt,freqChar] = cochlea_ieegequalpowerfilterbank ...
            (Fe, fMinHz, fMaxHz, nFilters, filterOrder, showResults)
% [bFilt,aFilt,freqChar] = COCHLEA_IEEGEQUALPOWERFILTERBANK 
%           (Fe, fMinHz, fMaxHz, nFilters, filtOrder, showResults)
% Computes the filter coefficients of the cochlea's filterbank using a
% frequency scale based on the mean intracerebral-EEG spectrum. This mean
% iEEG spectrum was estimated on a dozen of different iEEG files from
% different patients. A x->K/x.^b curve was fit. K and b values can be seen
% and changed in the iEEGequalpowerbands function. The frequency bands are
% then computed so that each band will have the same power on average. 
% The idea is to have the same number of output spikes for each channel to
% increase the STDP performance. 
% 
% The filterbank is an array of bandpass filter. The filters used are
% Butterworth filters. If any unstability is detected, the script gives a 
% warning the command line interface.
%
% INPUTS:
%   - Fe            : Sampling frequency (Hz)
%   - filtOrder     : Order of these filters (Butterworth filters)
%                    (must be multiple of 2)
%   - fMinHz        : Minimal frequency of analysis for the cochlea (Hz)
%   - fMaxHz        : Maximal frequency of analysis for the cochlea (Hz)
%   - nFilters      : Number of bandpass filters in the filterbank
%   - filtOrder     : Order of these filters (Butterworth filters)
%                    (must be multiple of 2)
%   - showResults   : If 1, show the frequency response of the filterbank
%
% OUTPUTS:
%   - bFilt         : Filter numerator coefficients
%   - aFilt         : Filter denominator coefficients
%   - freqChar      : Filter center frequencies (Hz)
%
% See also iEEGequalpowerbands, cochlea_eegfilterbank
%
% Author(s) : Martin Deudon (2016)

if nargin==5
    showResults = 0;
end

overlap     = 0.15;
fVect       = iEEGequalpowerbands(fMinHz,fMaxHz,nFilters);
freqBands   = [(1-overlap)*fVect(1:end-1),(1+overlap)*fVect(2:end)];

bFilt   = zeros(filterOrder+1,nFilters);
aFilt   = zeros(filterOrder+1,nFilters);
freqChar = zeros(1,nFilters);
for i=1:nFilters
    [bFilt(:,i),aFilt(:,i)] = butter(filterOrder/2,(2/Fe).*freqBands(i,:),'bandpass');
%     [b(:,i),a(:,i)] = ellip(filterOrder/2,1,80,(2/Fe).*freqBands(i,:));
    freqChar(i) = freqBands(i,1)+0.5*diff(freqBands(i,:));
    %- Check stability
    if ~isstable(bFilt(:,i),aFilt(:,i))
        warning(['Unstable filter detected for index ',num2str(i)]);
    end
end


%% Results
if showResults
    filterInd = 1:nFilters;
    figure;
    plot(filterInd,freqChar,'x'); hold on;
    fBandStart  = freqBands(:,1);
    fBandEnd    = freqBands(:,2);
    plot(repmat(filterInd,2,1),[fBandStart,fBandEnd]','color','r');
    xlabel('Filter Index'); ylabel('Caracteristic Frequency (Hz)');
    title('Caracteristic frequency of bandpass filters (CCIs) (iEEG equal power)');

    % Filter impulse responses
    nPoints = 50000;
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
    axis([min(freqBands(:))-5,max(freqBands(:))*2,-60,5]);
    xlabel('Frequency (Hz)'); ylabel('Gain (dB)'); 
    title(['Filter Bank Impulse Response - Filter order: ',num2str(filterOrder)]);
    ax(2) = subplot(212);
    semilogx(F,HdB);
    axis([min(freqBands(:)),Fe/2,-60,5]);
    xlabel('Frequency (Hz)'); ylabel('Gain (dB)'); 
    title(['Filter Bank Impulse Response - Filter order: ',num2str(filterOrder)]);
    linkaxes(ax,'x');
end


end

