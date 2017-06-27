
%% Parameters
eegFilepath     = 'C:\Users\deudon\Desktop\SAB\data\p25\SAB_Bipolaire\SAB_ENC_2.edf';
channelSelPos   = 94;
patientNumStr   = regexp(eegFilepath,[filesep,'p\d+'],'match');
patientNumStr   = patientNumStr{:};

%% Load signal
EEG     = pop_biosig(eegFilepath);
Fe      = EEG.srate;

%- Channel sel
for i=1:EEG.nbchan
    disp([num2str(i),' ',EEG.chanlocs(i).labels]);
end
x   = EEG.data(channelSelPos,:)';
disp(['Keeping channel(s) : ',EEG.chanlocs(channelSelPos).labels]);

%% Spectral Analysis (using Welch's power spectral density estimate)
nWindow     = round(length(x)/8);
overlap     = 0.5;
nOverlap    = round(nWindow*overlap);
nfft        = 16384;
fMinHz      = 3;    % Hz
fMaxHz      = 60;  % Hz

% [Pxx,f]     = pwelch(x,nWindow,nOverlap,nfft,Fe);
% [Pxx2,f2]   = periodogram(x,hanning(length(x)),nfft,Fe);
% figure;
% plot(f,(Pxx)); hold on;
% plot(f,500*(1./(f.^2)),'r');
% % plot(f2,10*log10(Pxx2),'r');
% xlim([fMinHz,fMaxHz]);
% xlabel('frequency (Hz)'); ylabel('Magnitude (dB)');


%% Mean Spectrum over all channels
eegChannelInd   = cellfun(@(x)strcmp(x(1:3),'EEG'),{EEG.chanlocs.labels});
nEegChannels    = sum(eegChannelInd);
EEGdata         = EEG.data(eegChannelInd,:);
PxxAll          = zeros(nEegChannels,1+nfft/2);
PxxAll2         = zeros(nEegChannels,1+nfft/2);
PxxAllconf      = zeros(nEegChannels,1+nfft/2,2);
for i=1:nEegChannels
    [PxxAll(i,:),f,PxxAllconf(i,:,:)] = pwelch(EEGdata(i,:),nWindow,nOverlap,nfft,Fe,'confidenceLevel',0.95);  
    PxxAll2(i,:)    = periodogram(EEGdata(i,:),hanning(length(x)),nfft,Fe);
end

PxxMean     = mean(PxxAll,1);
PxxConfMean = squeeze(mean(PxxAllconf));
PxxStd      = std(PxxAll);
PxxMean2    = mean(PxxAll2,1);
%- Curve fitting
fMinInd     = round((nfft/Fe)*fMinHz+1);
fMaxInd     = round((nfft/Fe)*fMaxHz+1);
PxxMeanSel  = PxxMean(fMinInd:fMaxInd);
PxxMeanSel2 = PxxMean2(fMinInd:fMaxInd);
fSel        = f(fMinInd:fMaxInd);
fitRes      = fit(fSel(:),PxxMeanSel(:),'power1');


%% Result

figure; hold on;
% plot(f,(PxxMean2),'r');
plot(f,(PxxMean));
fill([f(:);flipud(f(:))],[PxxMean(:)+PxxConfMean(:,1);flipud(PxxMean(:)-PxxConfMean(:,2))],'c','EdgeColor', [0.7,0.7,1],'FaceAlpha', 0.2);
%- Fitted curve
plot(f,fitRes.a*f.^fitRes.b,'k');
xlim([fMinHz,fMaxHz]);
xlabel('frequency (Hz)'); ylabel('Magnitude');
legend({'Mean PSD (PWelch)','95% CIs',sprintf('Fitted function: %.2f *f^%.2f',fitRes.a,fitRes.b)});
tempSep     = regexp(eegFilepath,filesep);
eegFilename = eegFilepath(tempSep(end)+1:end);
title([patientNumStr,' - Mean Power Spectral Density and Fitted function - signal : ',eegFilename]);

saveas(gcf,fullfile('meanSpectrum',[patientNumStr,'_',eegFilename(1:end-3),'.png']),'png');

%%
figure; hold on;
% plot(f,10*log10(PxxMean2),'r');
plot(f,10*log10(PxxMean));
%- Fitted curve
plot(f,10*log10(fitRes.a*f.^fitRes.b),'k');
fill([f(:);flipud(f(:))],[10*log10(PxxMean(:)+(PxxConfMean(:,1)));flipud(10*log10(PxxMean(:)-PxxConfMean(:,2)))],'c','EdgeColor', [0.7,0.7,1],'FaceAlpha', 0.2);
xlim([fMinHz,fMaxHz]);
xlabel('frequency (Hz)'); ylabel('Magnitude (dB)');
legend({'Mean PSD (PWelch)',sprintf('Fitted function: %.2f *f^%.2f',fitRes.a,fitRes.b)});
tempSep     = regexp(eegFilepath,filesep);
eegFilename = eegFilepath(tempSep(end)+1:end);
title([patientNumStr,' - Mean Power Spectral Density and Fitted function - signal : ',eegFilename]);

saveas(gcf,fullfile('meanSpectrum',[patientNumStr,'_',eegFilename(1:end-3),'_dB.png']),'png');

