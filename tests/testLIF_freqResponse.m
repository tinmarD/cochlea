% Frequency response of LIF neuron

Fe  = 2048; 
R   = 1;
C   = 0.1;

w = 0:0.01:pi;

ZMod = R./sqrt(1+(R*C*w).^2);


%% Results
figure;
plot(w,20*log10(ZMod));
