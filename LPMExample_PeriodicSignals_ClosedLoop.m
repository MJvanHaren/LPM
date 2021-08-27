clear all; 
close all; clc;
SetPlotLatexStyle;
rng(100);
[c1, c2, c3, c4, c5,c6,c7] = MatlabDefaultPlotColors();
%% system
s = tf('s');
zeta = 0.05;
Omega1 = 5; % rad/s
Ts = 1/250;
Omega2 = 4*Omega1; % rad/s
P = 5.6e3*(0.9*s^2+2*s+5.4e4)/(1.1*0.9*s^4+(1.1+0.9)*2*s^3+(1.1+0.9)*5.4e4*s^2);
C = 5*(1/(10*2*pi)*s+1)/(1/(40*2*pi)*s+1);
% P = [P 0;0 P];
% C = [C 0;0 C];
P = c2d(P,Ts,'tustin');
C = c2d(C,Ts,'tustin');
Nu = size(C,2);
Ny = size(P,1);

n = 4; % window size left and right
R=2; % max. degree of polynomial in LPM
tperiod = 10; % time for signal
Per = 20; % amount of periods
%% signal
Np = tperiod/Ts; % amount of samples in period
Ns = Per*Np; % total amount of samples
N = floor(Np/2);
f = linspace(0, 1 - 1/N, N) * (1/Ts)/2; % available frequencies:
A = ones(N,1); % amplitude distribution
tp = (0:Ts:(Np-1)*Ts)';
t = (0:Ts:(Ns-1)*Ts)';

% custom multisine
rp = zeros(Np,1); % periodic signal
rp2 = zeros(Np,1); % periodic signal
for k = 1:N
   rp = rp+A(k)*sin(2*pi*f(k)*tp+rand*2*pi);
   rp2 = rp2+A(k)*sin(2*pi*f(k)*tp+rand*2*pi);
end
r = repmat(rp,Per,1);
% r = [r repmat(rp2,Per,1)];

% white noise
% r = randn(Ns,Nu);
%% simulate and LPM/ETFE
y = lsim(P*C/(1+P*C),r,t); % Pintelon 2012 Figure 7-4!
u = lsim(C/(1+P*C),r,t); % Pintelon 2012 Figure 7-4!

[P_LPM] = LPMClosedLoopPeriodicRobustFRM(u,y,r,n,R,Per);
%% plotting
figure(1); clf;
magP = bode(P,f*2*pi); hold on;
for i = 1:Nu
    for ii = 1:Ny
        subplot(Ny,Nu,i+(ii-1)*Nu);
        semilogx(f,20*log10(abs(squeeze(magP(ii,i,:)))),'Color',c1); hold on;
        semilogx(f,20*log10(abs(squeeze(P_LPM(ii,i,:)))),'Color',c2);
        set(gca,'xscale','log');
        xlim([f(2) f(end)]);
    end
end
% semilogx(f,20*log10(abs(T_LPM(1:length(f)))),'Color',c3);
legend('True plant','Estimated plant')

% figure(2); clf;
% [magG,~,~] = bode(P,f*2*pi);
% absG = (squeeze(magG));
% semilogx(f,mag2db(abs(abs(P_LPM')-absG)),'Color',c2); hold on 
% % [magG,~,~] = bode(P,P_ETFE.Frequency*128*2);
% % absG = (squeeze(magG));
% % semilogx(P_ETFE.Frequency*128/pi,mag2db(abs(squeeze(abs(P_ETFE.ResponseData))-absG)),'Color',c4);
% xlabel('Frequency [Hz]'); xlim([f(1) f(end)]);
% ylabel('Estimation Error [dB]');
% legend('LPM','ETFE','location','best')