% (i)EEG mean spectrum estimation - Reads all the edf files in the input
% folder, estimated the power spectrum using Welch's method. Compute the
% mean of all the power spectrum estimates. Fit a power law (x->K*x^b) 
% curve on this data. Curve fitting is done between fMinHz and fMaxHz (Hz)
%
% Author(s) : Martin Deudon (2016)

%% Parameters
dataDirPath = 'C:\Users\deudon\Desktop\SAB\_Data\COG_015\SAB_Bipolaire';
resDirPath  = 'C:\Users\deudon\Desktop\M4\_Results\meanSpectrum';
nfft        = 16384;    % Number of point in the fft
fMinHz      = 5;        % Curve fitting is done between fMinHz and fMaxHz (Hz)
fMaxHz      = 45;       % Curve fitting is done between fMinHz and fMaxHz (Hz)

%% 
fileList    = rdir(fullfile(dataDirPath,'**',filesep,'*.edf'));
nFiles      = length(fileList);
disp([num2str(nFiles),' EDF files found']);

%% 
PxxMeanAll      = zeros(nFiles,1+nfft/2);
PxxMeanConfAll  = zeros(nFiles,1+nfft/2,2);
fHzToInd        = @(fHz,Fe,nfft)round(1+nfft*fHz/Fe);
if nFiles==0
    disp(['Could not find any .EDF files in ',dataDirPath]);
    return;
end
for iFile=1:nFiles
	EEG  = pop_biosig(fileList(iFile).name);
%     if  EEG.srate>256
%         disp('Resampling...');
%         EEG = pop_resample(EEG,256);
%     end
    Fe = EEG.srate;
    nWindow     = round(EEG.pnts/8);
    overlap     = 0.5;
    nOverlap    = round(nWindow*overlap);
    eegChannelInd   = cellfun(@(x)strcmp(x(1:3),'EEG'),{EEG.chanlocs.labels});
    nEegChannels    = sum(eegChannelInd);
    EEGdata         = EEG.data(eegChannelInd,:);
    PxxAll          = zeros(nEegChannels,1+nfft/2);
    PxxAllconf      = zeros(nEegChannels,1+nfft/2,2);
    for iChan=1:nEegChannels
        [PxxAll(iChan,:),f,PxxAllconf(iChan,:,:)] = pwelch(EEGdata(iChan,:),nWindow,nOverlap,nfft,Fe,'confidenceLevel',0.95);  
    end
    PxxMeanAll(iFile,:)         = mean(PxxAll,1);
    PxxMeanConfAll(iFile,:,:)   = squeeze(mean(PxxAllconf));
    disp([num2str(iFile),'/',num2str(nFiles)]);
end
PxxMean     = mean(PxxMeanAll);
PxxConfMean = squeeze(mean(PxxMeanConfAll));

%% Curve fitting
fMinInd     = round((nfft/Fe)*fMinHz+1);
fMaxInd     = round((nfft/Fe)*fMaxHz+1);
PxxMeanSel  = PxxMean(fMinInd:fMaxInd);
fSel        = f(fMinInd:fMaxInd);
fitRes      = fit(fSel(:),PxxMeanSel(:),'power1');


%% Results
figure; hold on;
% plot(f,(PxxMean2),'r');
plot(f,(PxxMean));
fill([f(:);flipud(f(:))],[PxxMean(:)+PxxConfMean(:,1);flipud(PxxMean(:)-PxxConfMean(:,2))],'c','EdgeColor', [0.7,0.7,1],'FaceAlpha', 0.2);
%- Fitted curve
plot(f,fitRes.a*f.^fitRes.b,'k');
xlim([fMinHz,fMaxHz]);
xlabel('frequency (Hz)'); ylabel('Magnitude');
legend({'Mean PSD (PWelch)','95% CIs',sprintf('Fitted function: %.2f *f^%.2f',fitRes.a,fitRes.b)});
title('P27 - Mean Power Spectral Density and Fitted function');

% saveas(gcf,fullfile(resDirPath,'P27_meanSpectrumAcross.png'),'png');
