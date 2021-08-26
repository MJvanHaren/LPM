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
G0s = P/(1+P*C);
G0 = c2d(G0s,Ts,'tustin');
ssG0 = ss(G0);

n = 6 ; % window size left and right
R=2;

%% signal
tperiod = 40;
Np = tperiod/Ts;
N = floor(Np/2);
T0 = Ts*Np;
% available frequencies:
k = (1:N)'; 
f = linspace(0, 1 - 1/N, N) * (1/Ts)/2;
A = ones(N,1); % amplitude distribution
t = (0:Ts:(Np-1)*Ts)';

% custom multisine
u = zeros(Np,1);
for k = 1:N
   u = u+A(k)*sin(2*pi*f(k)*t+rand*2*pi);
end

% white noise
% u = randn(Np,1);

figure;
plot(t,u);


%% simulate and LPM/ETFE
% y = lsim(G0,u,t)+0.001*randn(Np,1);
y = lsim(G0,u,t);
G_ETFE = etfe([y u],100,N);
[G_LPM,T_LPM] = LPMOpenLoop(u,y,n,R);
%%
figure
opts = bodeoptions;
opts.FreqUnits = 'Hz';
opts.xlim = [1e-1 125];
bodemag(G0,opts); hold on;
semilogx(f,20*log10(abs(G_LPM)),'Color',c2); hold on %-(n-1)*Ts or -n*Ts is perfect! -> fix in freq window?
semilogx(f,20*log10(abs(T_LPM(1:length(f)))),'Color',c3);
semilogx(G_ETFE.Frequency*128/pi,20*log10(squeeze(abs(G_ETFE.ResponseData))),'Color',c4);
legend('True plant','Estimated plant','Estimated transient plant','ETFE')

figure
[magG,~,~] = bode(G0,f*2*pi);
absG = (squeeze(magG));
semilogx(f,mag2db(abs(abs(G_LPM)-absG)),'Color',c2); hold on 
[magG,~,~] = bode(G0,G_ETFE.Frequency*128*2);
absG = (squeeze(magG));
semilogx(G_ETFE.Frequency*128/pi,mag2db(abs(squeeze(abs(G_ETFE.ResponseData))-absG)),'Color',c4);
xlabel('Frequency [Hz]'); xlim([f(1) f(end)]);
ylabel('Estimation Error [dB]');
legend('LPM','ETFE','location','best')



