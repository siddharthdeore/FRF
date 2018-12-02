%% -------------- Frequency Response Function estimate -------------- %%
%%
%% This example evaluates the FRF of a system with two degrees of freedom
%% excited by white noise. Additive uncorrelated measurement noise is
%% considered acting on:
%% case a) the output only, 
%% case b) the input only, 
%% case c) both the input and the output. 
%% Three estimators, namely H1 H2 and Hv, are compared for each case. 
%%
% boot commands:
clear all; close all;clc 
% Program:
% Define the system parameters:
A1 = 400; 
A2 = 60; 
f1 = 25; 
f2 = 35; 
wn1 = 2*pi*f1; 
wn2 = 2*pi*f2; 
zeta1 = 0.001; 
zeta2 = 0.003;
wd1 = sqrt(1-zeta1^2)*wn1; 
wd2 = sqrt(1-zeta2^2)*wn2;
% Define the acquisition parameters:
fs = 100; 
T1 = 10; 
t1 = 0:1/fs:T1-1/fs;
% Evaluate the system impulse response function h(t) and the FRF H(f)
h = (A1/wd1)*exp(-zeta1*wn1*t1).*sin(wd1*t1)+(A2/wd2)*exp(-zeta2*wn2*t1).*sin(wd2*t1);
H = fft(h);
% Generate the input signal x(t), as white noise:
T = 100;% -> 100 segments without overlap
x = randn(1,T*fs); 
% Generate the output signal y(t) as the system response to the input
y = filter(h,1,x); % we do not scale for convenience.
% Generate the measuremet noise acting on input and output for case (c)
nx = 0.5*randn(size(x)); % nx=0 for case (a)
ny = 0.5*randn(size(y)); % ny=0 for case (b)
% Calculate the (one-sided) spectral density functions using the segment
% averaging method with hanning window, 50% overlap and 1000 segments. 
% Then calculate the FRF and the coherence function.  
disp('Wait patiently...');
% Case a)
ya = y+ny;
[Gxx, f] = cpsd(x(1:T*fs),x(1:T*fs), hanning(T1*fs),T1*fs/2, T1*fs, fs);
[Gyy, f] = cpsd(ya(1:T*fs),ya(1:T*fs), hanning(T1*fs),T1*fs/2, T1*fs, fs);
[Gxy, f] = cpsd(x(1:T*fs),ya(1:T*fs), hanning(T1*fs),T1*fs/2, T1*fs, fs);
[Gyx, f] = cpsd(ya(1:T*fs),x(1:T*fs), hanning(T1*fs),T1*fs/2, T1*fs, fs);
H1a = Gxy./Gxx; 
H2a = Gyy./Gyx;
HTa = (Gyy-Gxx + sqrt((Gxx-Gyy).^2 + 4*abs(Gxy).^2))./(2*Gyx);
Gammaa = abs(Gxy).^2./(Gxx.*Gyy);
disp('                ...still waiting?');
% case b)
xb = x+nx; 
[Gxx, f] = cpsd(xb(1:T*fs),xb(1:T*fs), hanning(T1*fs),T1*fs/2, T1*fs, fs);
[Gyy, f] = cpsd(y(1:T*fs),y(1:T*fs), hanning(T1*fs),T1*fs/2, T1*fs, fs);
[Gxy, f] = cpsd(xb(1:T*fs),y(1:T*fs), hanning(T1*fs),T1*fs/2, T1*fs, fs);
[Gyx, f] = cpsd(y(1:T*fs),xb(1:T*fs), hanning(T1*fs),T1*fs/2, T1*fs, fs);
H1b = Gxy./Gxx; 
H2b = Gyy./Gyx;
HTb = (Gyy-Gxx + sqrt((Gxx-Gyy).^2 + 4*abs(Gxy).^2))./(2*Gyx);
Gammab = abs(Gxy).^2./(Gxx.*Gyy);
disp('                                  Time for a coffee?');
disp('                                                        S')
disp('                                                       __s__')
disp('                                                      |     |>')
disp('                                                      \_____/')
% case c)
% xc=xb; yc=ya;
[Gxx, f] = cpsd(xb(1:T*fs),xb(1:T*fs), hanning(T1*fs),T1*fs/2, T1*fs, fs);
[Gyy, f] = cpsd(ya(1:T*fs),ya(1:T*fs), hanning(T1*fs),T1*fs/2, T1*fs, fs);
[Gxy, f] = cpsd(xb(1:T*fs),y(1:T*fs), hanning(T1*fs),T1*fs/2, T1*fs, fs);
[Gyx, f] = cpsd(ya(1:T*fs),xb(1:T*fs), hanning(T1*fs),T1*fs/2, T1*fs, fs);
H1c = Gxy./Gxx; 
H2c = Gyy./Gyx;
HTc = (Gyy-Gxx + sqrt((Gxx-Gyy).^2 + 4*abs(Gxy).^2))./(2*Gyx);
Gammac = abs(Gxy).^2./(Gxx.*Gyy);
clc;disp('Espresso');
% Plot of the results
figure
subplot(2,1,1),plot(f,20*log10(abs(H1a)), f,20*log10(abs(H2a)),f,20*log10(abs(HTa)),f,20*log10(abs(H(1:length(f)))),':')
title('case a: output noise only')
legend('|\itH\rm_1|','|\itH\rm_2|','|\itH\rm_T|','|\itH\rm|')
xlabel('Frequency (Hz)'); ylabel('Modulus (dB)')
axis([0 25 -35 25])
subplot(2,1,2),plot(f,-unwrap(angle(H1a)), f,-unwrap(angle(H2a)),f,-unwrap(angle(HTa)),f,unwrap(angle(H(1:length(f)))),':')
xlabel('Frequency (Hz)'); ylabel('Phase (rad)')
axis([0 25 -4 0.5])
legend('\angle H_1','\angle H_2','\angle H_T','\angle H')
%
figure
subplot(2,1,1),plot(f,20*log10(abs(H1b)), f,20*log10(abs(H2b)),f,20*log10(abs(HTb)),f,20*log10(abs(H(1:length(f)))),':')
title('case b: input noise only')
legend('|\itH\rm_1|','|\itH\rm_2|','|\itH\rm_T|','|\itH\rm|')
xlabel('Frequency (Hz)'); ylabel('Modulus (dB)')
axis([0 25 -35 25])
subplot(2,1,2),plot(f,-unwrap(angle(H1b)), f,-unwrap(angle(H2b)),f,-unwrap(angle(HTb)),f,unwrap(angle(H(1:length(f)))),':')
xlabel('Frequency (Hz)'); ylabel('Phase (rad)')
axis([0 25 -4 0.5])
legend('\angle H_1','\angle H_2','\angle H_T','\angle H')
%
figure
subplot(2,1,1),plot(f,20*log10(abs(H1c)), f,20*log10(abs(H2c)),f,20*log10(abs(HTc)),f,20*log10(abs(H(1:length(f)))),':')
title('case c: both input and output noise')
legend('|\itH\rm_1|','|\itH\rm_2|','|\itH\rm_T|','|\itH\rm|')
xlabel('Frequency (Hz)'); ylabel('Modulus (dB)')
axis([0 25 -35 25])
subplot(2,1,2),plot(f,-unwrap(angle(H1c)), f,-unwrap(angle(H2c)),f,-unwrap(angle(HTc)),f,unwrap(angle(H(1:length(f)))),':')
xlabel('Frequency (Hz)'); ylabel('Phase (rad)')
axis([0 25 -4 0.5])
legend('\angle H_1','\angle H_2','\angle H_T','\angle H')
%
figure
plot(f, Gammaa,f,Gammab,f,Gammac); 
xlabel('Frequency (Hz)'); ylabel('Coherence function')
legend('case a','case b','case c')
axis([0 25 0 1])
clc
% END
