% Filter design similar to O. Bichlet thesis

%% Parameters
Fe          = 41000;
fMinHz      = 20;
fMaxHz      = 20000;
nFilter     = 300;
filtOrder   = 4;

Qear        = 9.26449;
BWmin       = 24.7;

%% 
tic
%- Calcul of the central frequencies of band-pass filters
filterFcHz      = erbspace(fMinHz,fMaxHz,nFilter);
ERB             = 24.7+filterFcHz/9.265;

%- Filter design
% Hypothesis of butterworth band-pass filters with pass-band centered on
% the caracteristic frequency
fBandStart  = max(fMinHz,filterFcHz-ERB/2);
fBandEnd    = min(fMaxHz,filterFcHz+ERB/2);
fBandStart(fBandStart<0)  = 0;
fBandEnd  (fBandEnd>Fe/2) = Fe/2.01;
aFilt       = zeros(nFilter,filtOrder+1);
bFilt       = zeros(nFilter,filtOrder+1);
for i=1:nFilter
    [bFilt(i,:),aFilt(i,:)] = butter(filtOrder/2,(2/Fe)*[fBandStart(i),fBandEnd(i)],'bandpass');
end
toc;

%% Results
filterInd   = 1:nFilter;

figure;
plot(filterInd,filterFcHz,'x'); hold on;
plot(repmat(filterInd,2,1),[fBandStart;fBandEnd],'color','r');
xlabel('Filter Index'); ylabel('Caracteristic Frequency (Hz)');
title('Caracteristic frequency of bandpass filters (CCIs) with ERB');

% Filter impulse responses
nPoints = 1000;
H       = zeros(nFilter,nPoints);
F       = zeros(1,nPoints);
s.xunits = 'Hz';
s.yunits = 'db';
s.plot   = 'mag';
s.fvflag = 0;
s.yphase = 'degrees';
for i=1:nFilter
    [H(i,:),F] = freqz(bFilt(i,:),aFilt(i,:),nPoints,Fe);
end
HdB = 20*log10(abs(H));

figure;
plot(F,HdB);
axis([fMinHz,Fe/2,-60,5]);
xlabel('Frequency (Hz)'); ylabel('Gain (dB)'); 
title(['Filter Bank Impulse Response - Filter order: ',num2str(filtOrder)]);

% figure;
% semilogx(F,HdB);
% axis([fMinHz,Fe/2,-60,5]);
% xlabel('Frequency (Hz)'); ylabel('Gain (dB)'); 
% title(['Filter Bank Impulse Response - Filter order: ',num2str(filtOrder)]);

