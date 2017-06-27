clear all;

%% Parameters
wavFilepath     = 'C:\Users\deudon\Desktop\M4\_Data\sounds\cigars.wav';
% wavFilepath     = 'C:\Users\deudon\Desktop\M4\_Data\sounds\testSounds\square440hz.wav';
% wavFilepath     = 'C:\Users\deudon\Desktop\M4\_Data\sounds\testSounds\R10N1.wav';
resultsFigDir   = 'C:\Users\deudon\Desktop\M4\Resultats\Res_main_LIFmodel';
saveFigure      = 0;

fMinHz          = 20;
fMaxHz          = 8000;
nFilters        = 300;
%- rectifier type
rectifierType   = 'classic_v2';    % 'classic' or 'hybrid'
% rectifierType   = 'hybrid';     % 'classic' or 'hybrid'
%- LP filter
filterOrderLP   = 2;
freqCharRatio   = 1.5;
%- filter type
bankFilterType  = 'linear'; % 'linear' or 'erb'
% bankFilterType  = 'erb'; % 'linear' or 'erb'
filterOrder 	= 4;    % Only multiple of 2 order
simDuration     = 1;

seppos          = regexp(wavFilepath,filesep);
wavFilename     = wavFilepath(seppos(end)+1:end-4);


%% Open wav file
[x,Fe]          = audioread(wavFilepath,'double');
if size(x,2)~=1;
    x = x(:,1);
end
if length(x)>simDuration*Fe
    x = x(1:simDuration*Fe);
end
if Fe<=2*fMaxHz; error('fMaxHz must be inferior to half the sampling frequency'); end;
xDurationSec    = length(x)/Fe;
disp(['Sampling frequency: ',num2str(Fe),' Hz.']);
disp(['Signal duration: ',num2str(xDurationSec),' sec.']);
disp(['min(x) : ',num2str(min(x)),char(10),'Max(x) : ',num2str(max(x))]);


%% Normalization of audio file
a       = 2/(max(x)-min(x));
b       = 1-a*max(x);
xNorm   = a*x+b;
% x       = xNorm;
% disp('Normalizing to [-1,1]');


%% Filter Design
tic;
if strcmpi(bankFilterType,'erb')
    [b,a,freqChar,ERB]  = cochlea_getfilterbank(Fe,fMinHz,fMaxHz,nFilters,filterOrder);
elseif strcmpi(bankFilterType,'linear')
    [b,a,freqChar]      = cochlea_linearfilterbank(Fe,fMinHz,fMaxHz,nFilters,filterOrder);
else
    error('You must specify bankFilterType (erb or linear)');
end
filterDesignTime    = toc;


%% Bank filtering with rectifier
if strcmp(rectifierType,'classic')
    [bLPfilt,aLPfilt]   = butter(1,2*65/Fe,'low');
    [Yfilt,Yrect]       = cochlea_bankfilterrectifier_classic(b,a,bLPfilt,aLPfilt,x);
elseif strcmp(rectifierType,'classic_v2')   
        aLPfilt = zeros(filterOrderLP+1,nFilters);
        bLPfilt = zeros(filterOrderLP+1,nFilters);
        for iFilt=1:nFilters
            [bLPfilt(:,iFilt), aLPfilt(:,iFilt)] = butter(filterOrderLP,(2/Fe)*freqChar(iFilt)/freqCharRatio);
        end      
        [Yfilt,Yrect]       = cochlea_bankfilterrectifier_classic_v2(b,a,bLPfilt,aLPfilt,x);
elseif strcmp(rectifierType,'hybrid')
    %- Characteristic Period (for the rectifier)
    periodChar          = round(Fe./freqChar);
    periodChar          = periodChar(:);
    [Yfilt,Yrect]       = cochlea_bankfilterrectifier_hybrid(b,a,x,periodChar);
else
    error('You must specify rectifierType (erb or linear)');
end

%- Compression
Ycomp = Yrect.^(1/3);


%% Visualisation
t = linspace(0,xDurationSec,length(x));
tsec2ind  	= @(x)(x-t(1))*Fe+1;
xlims = [0.1,0.13];
ylims = [-0.1,0.1];
iChan = 42;

%- Hybrid rectifier
% figure;
% ax(1) = subplot(311);
% plot(t,x); title(['Raw input signal : ',wavFilename]);
% ylim([-1,1]);
% ax(2) = subplot(312);
% plot(t,Yfilt(:,iChan)); title(['Filtered signal - channel ',num2str(iChan),' : ',num2str(round(freqChar(iChan))),' Hz']);
% ylim(ylims);
% ax(3) = subplot(313);
% plot(t,Yrect_hyb(:,iChan)); 
% ylim(ylims);
% title(['Filtered/rectified signal (Hybrid Rectification) - channel ',num2str(iChan),' : ',num2str(round(freqChar(iChan))),' Hz']);
% linkaxes(ax,'x');
% xlim(xlims);
% % 
%- Classic rectifier
figure;
ax(1) = subplot(311);
plot(t,x); title(['Raw input signal : ',wavFilename]);
ylim([-1,1]);
ax(2) = subplot(312);
plot(t,Yfilt(:,iChan)); title(['Filtered signal - channel ',num2str(iChan),' : ',num2str(round(freqChar(iChan))),' Hz']);
hold on;
Yfilt_hwrect = Yfilt(:,iChan); Yfilt_hwrect(Yfilt_hwrect<0)=0;
plot(t,Yfilt_hwrect,'r:');
plot(t,Yrect(:,iChan),'g'); 
ylim(ylims);
legend({'Filtered signal','Half-wave rectified','Rectified signal'});
ax(3) = subplot(313);
plot(t,Ycomp(:,iChan)); 
% ylim(ylims);
title(['Half-wave rectified + Low-Pass filter (Classic Rectification) - channel ',num2str(iChan),' : ',num2str(round(freqChar(iChan))),' Hz']);
linkaxes(ax,'x');
axis tight;




%% - Filtered and Rectified signal
t = linspace(0,xDurationSec,length(x));
tsec2ind  	= @(x)(x-t(1))*Fe+1;

figure;
imagesc(Yrect');
colorbar;
set(gca,'Ydir','normal');
tStepSec    = 0.2;
tTickVal 	= 0:tStepSec:floor(xDurationSec);
tTickInd    = tsec2ind(tTickVal);
set(gca,'xtick',tTickInd,'xticklabel',num2str(tTickVal(:)));
ytick       = get(gca,'ytick');
yTickVal    = round(freqChar(ytick));
set(gca,'ytick',ytick,'yticklabel',num2str(yTickVal(:)));
xlabel('time (s)'); ylabel('characteristic frequency (Hz)');
title(['Rectified signal (',rectifierType,' Rectification)']);
if saveFigure
    saveas(gcf,fullfile(resultsFigDir,[wavFilename,'_Filt.png'])); 
end


%% - Compressed signal
% 
figure;
imagesc(Ycomp');
colorbar;
set(gca,'Ydir','normal');
tStepSec    = 0.2;
tTickVal 	= 0:tStepSec:floor(xDurationSec);
tTickInd    = tsec2ind(tTickVal);
set(gca,'xtick',tTickInd,'xticklabel',num2str(tTickVal(:)));
ytick       = get(gca,'ytick');
yTickVal    = round(freqChar(ytick));
set(gca,'ytick',ytick,'yticklabel',num2str(yTickVal(:)));
xlabel('time (s)'); ylabel('characteristic frequency (Hz)');
title(['Compressed signal (',rectifierType,' Rectification)']);
if saveFigure
    saveas(gcf,fullfile(resultsFigDir,[wavFilename,'_Filt.png'])); 
end
