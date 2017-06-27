function [] = testFiltrage(wavFilepath, resultsFigDir, saveFigure)


%% Parameters
if nargin==0
    % wavFilename     = 'C:\Users\deudon\Desktop\M4\Prog\mexfile\cigars.wav';
    wavFilepath     = 'C:\Users\deudon\Desktop\M4\sounds\testSounds\whitenoise.wav';
    % wavFilename     = 'C:\Users\deudon\Desktop\M4\sounds\Gaussian_sounds\10repetitions\R10N1.wav';
end
if nargin<2
    resultsFigDir   = 'C:\Users\deudon\Desktop\M4\Resultats\Res_main_LIFmodel';
end
if nargin<3
    saveFigure      = 0;
elseif saveFigure==1
    %- External call with figure saving
end
fMinHz          = 20;
fMaxHz          = 20000;
nFilters        = 1024;
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
[b,a,freqChar]  = cochlea_getfilterbank(Fe,fMinHz,fMaxHz,nFilters,filterOrder);
%- Characteristic Period (for the rectifier)
periodChar      = round(Fe./freqChar);
periodChar      = periodChar(:);


%% Bank filtering with rectifier
tic;
[Yfilt,Yrect] = cochlea_bankfilterrectifier(b,a,x,periodChar);
Yfiltcomp = Yfilt;
Yfiltcomp(Yfiltcomp<0) = 0;
Yfiltcomp = Yfiltcomp.^(1/3);

Yfiltrectcomp = Yrect;
Yfiltrectcomp(Yfiltrectcomp<0) = 0;
Yfiltrectcomp = Yfiltrectcomp.^(1/3);

[b,a]           = butter(1,2*65/Fe,'low');
Yfiltrect2      = Yfilt;
Yfiltrect2(Yfiltrect2<0) = 0;
Yfiltrect2      = filter(b,a,Yfiltrect2);
Yfiltrectcomp2  = Yfiltrect2.^(1/3);

%% - Filtered and Rectified signal
t = linspace(0,xDurationSec,length(x));
tsec2ind  	= @(x)(x-t(1))*Fe+1;

figure;
imagesc(Yfilt');
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
title(['Filtered signal - ',num2str(nFilters),' filters (order ',num2str(filterOrder),') - fmin: ',num2str(fMinHz),'Hz - fMax: ',num2str(fMaxHz),'Hz - FE: ',num2str(Fe),'Hz']);
if saveFigure
    saveas(gcf,fullfile(resultsFigDir,[wavFilename,'_Filt.png'])); 
end

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
title(['Filt/Rectified signal - ',num2str(nFilters),' filters (order ',num2str(filterOrder),') - fmin: ',num2str(fMinHz),'Hz - fMax: ',num2str(fMaxHz),'Hz - FE: ',num2str(Fe),'Hz']);
if saveFigure
    saveas(gcf,fullfile(resultsFigDir,[wavFilename,'_FiltRect.png'])); 
end

figure;
imagesc(Yfiltcomp');
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
title(['Filt/Compressed signal - ',num2str(nFilters),' filters (order ',num2str(filterOrder),') - fmin: ',num2str(fMinHz),'Hz - fMax: ',num2str(fMaxHz),'Hz - FE: ',num2str(Fe),'Hz']);
if saveFigure
    saveas(gcf,fullfile(resultsFigDir,[wavFilename,'_FiltComp.png'])); 
end

figure;
imagesc(Yfiltrectcomp');
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
title(['Filt/Rect/Compressed signal - Hybrid rectifier - ',num2str(nFilters),' filters (order ',num2str(filterOrder),') - fmin: ',num2str(fMinHz),'Hz - fMax: ',num2str(fMaxHz),'Hz - FE: ',num2str(Fe),'Hz']);
if saveFigure
    saveas(gcf,fullfile(resultsFigDir,[wavFilename,'_FiltRectComp_hybrid.png'])); 
end


figure;
imagesc(Yfiltrectcomp2');
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
title(['Filt/Rect/Compressed signal - low-pass-filter rectifier - ',num2str(nFilters),' filters (order ',num2str(filterOrder),') - fmin: ',num2str(fMinHz),'Hz - fMax: ',num2str(fMaxHz),'Hz - FE: ',num2str(Fe),'Hz']);
if saveFigure
    saveas(gcf,fullfile(resultsFigDir,[wavFilename,'_FiltRectComp_lowpass.png'])); 
end
end