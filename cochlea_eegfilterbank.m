function [bFilt,aFilt,freqChar]  = cochlea_eegfilterbank ...
    (Fe,filterOrder,showResults,freqBands)
% [bFilt,aFilt,freqChar]  = COCHLEA_EEGFILTERBANK 
%               (Fe,filterOrder, showResults, freqBands)
% Computes the filter coefficients of the cochlea's filterbank using the 
% common EEG bands. By default, these frequency bands are (Hz): 
%   [0.5,3;3,7;7,12;12,15;15,20;20,30;30,40;40,60;60,90]
% Other bands can be specified by setting the last argument freqBands as a
% 2D vector [nBands,2]
% The filterbank is an array of bandpass filter. The filters used are
% Butterworth filters. If any unstability is detected, the script gives a 
% warning the command line interface.
%
% INPUTS:
%   - Fe            : Sampling frequency (Hz)
%   - filtOrder     : Order of these filters (Butterworth filters)
%                    (must be multiple of 2)
%   - showResults   : If 1, show the frequency response of the filterbank
%   - freqBands     : Limits of the frequency bands (Hz) - 2D vector [nBands,2]
%
% OUTPUTS:
%   - bFilt         : Filter numerator coefficients
%   - aFilt         : Filter denominator coefficients
%   - freqChar      : Filter center frequencies (Hz)
%
% See also cochlea_ieegequalpowerfilterbank
%
% Author(s) : Martin Deudon (2016)

if nargin<3
    showResults = 0;
end
if nargin<4
    freqBands = [0.5,3;3,7;7,12;12,15;15,20;20,30;30,40;40,60;60,90];
end
nFilters  = size(freqBands,1);

bFilt   = zeros(filterOrder+1,nFilters);
aFilt   = zeros(filterOrder+1,nFilters);
freqChar = zeros(1,nFilters);
for i=1:nFilters
    [bFilt(:,i),aFilt(:,i)] = butter(filterOrder/2,(2/Fe).*freqBands(i,:),'bandpass');
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