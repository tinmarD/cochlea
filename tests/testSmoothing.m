
% Smoothing function (LP filter) of the form : y(n+1) = a*y(n)+(1-a)*x(n)
% with 0<a<1
% Determine the frequency response and the -3db cutoff frequency

%- Calcul of the filter gain given a and w (pulsation)
fMod    = @(w,a)abs(1-a)./sqrt(1-2*a*cos(w)+a.^2);
%- Calcul of the cutoff frequency given a and b (value of cutoff)
fCutoff = @(a,b)acos((b^2*(1+a^2)-(1-a)^2)/(2*b^2*a));
Fe      = 44100; % Hz


%% Results 1
a       = 0.15;
w       = 0:0.001:pi;
wc      = fCutoff(a,10^-(3/20));
fc      = Fe*wc/(2*pi);
f       = Fe*w/(2*pi);

figure;
plot(f,20*log10(fMod(w,a)),'b'); hold on;
axis tight;
plot(xlim,[-3,-3],'c-.');
plot([fc,fc],ylim,'c-.');
xlabel(['f (Hz) (Fe = ',num2str(Fe),'Hz)']); ylabel('Gain (dB)');
title(['Smoothing frequency response - a=',num2str(a)]);


%% Results 2 
aValues = 0:0.01:1;
yGain   = 10^-(3/20);
wc      = zeros(1,length(aValues));
for i=1:length(aValues)
    a       = aValues(i);
    wc(i)   = fCutoff(a,yGain);
    if ~isreal(wc(i)); wc(i)=0; end;
end
fc = Fe*wc/(2*pi);
figure;
plot(aValues,fc);
xlabel('smoothing parameter value (a)');
ylabel('-3db cutoff frequency (Hz)');
    



