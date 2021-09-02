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
Nu = size(P,2);
Ny = size(P,1);

n = 4; % window size left and right
R=2; % max. degree of polynomial in LPM
tperiod = 40; % time for signal
%% signal
Np = tperiod/Ts;
N = floor(Np/2);
f = linspace(0, 1 - 1/N, N) * (1/Ts)/2; % available frequencies:
A = ones(N,1); % amplitude distribution
t = (0:Ts:(Np-1)*Ts)';

% custom multisine
r = zeros(Np,1); % periodic signal
rp2 = zeros(Np,1); % periodic signal
for k = 1:N
   r = r+A(k)*sin(2*pi*f(k)*t+rand*2*pi);
   rp2 = rp2+A(k)*sin(2*pi*f(k)*t+rand*2*pi);
end
if Ny==2
    r = [r rp2];
end

% white noise
% r = randn(Np,1);
%% simulate and LPM/ETFE
L = P*C; Li = C*P;
y = lsim(L/(1+L),r,t); % Pintelon 2012 Figure 7-4!
u = lsim(C/(1+Li),r,t); % Pintelon 2012 Figure 7-4!
% G_ETFE = zeros(Ny,Nu,N);
% P_ETFE = zeros(Ny,Nu,N);
% for i = 1:Nu
%     for ii = 1:Ny
%         G_ETFE(ii,i,:) = etfe([y(:,ii) u(:,i)],50,N).ResponseData;
%         
%     end
% end
% [G_ETFE,f_ETFE] = modalfrf(u,y,1/Ts,500);
% [magC,~,~] = bode(C,f_ETFE*2*pi);
% for k = 1:length(f_ETFE)
%     invC = inv(magC(:,:,k));
%     P_ETFE(:,:,k) = (inv(invC*reshape(G_ETFE(k,:,:),Ny,Nu)')-eye(Nu))*invC;
% end
% f_ETFE = etfe([y(:,1) u(:,1)],50,N).Frequency*128/pi;
[P_LPM] = LPMClosedLoopArbitrary(u,y,r,n,R);
%% plotting
figure(1); clf;
magP = bode(P,f*2*pi);
for i = 1:Nu
    for ii = 1:Ny
        subplot(Ny,Nu,i+(ii-1)*Nu);
        semilogx(f,20*log10(abs(squeeze(magP(ii,i,:)))),'Color',c1); hold on;
        semilogx(f,20*log10(abs(squeeze(P_LPM(ii,i,:)))),'Color',c2);
%         semilogx(f_ETFE,20*log10(abs(squeeze(G_ETFE(:,ii,i)))),'Color',c4);
        set(gca,'xscale','log');
        xlim([f(2) f(end)]);
    end
end
legend('True plant','Estimated plant','ETFE')

figure(2); clf;
difP = abs(abs(P_LPM)-magP);
for i = 1:Nu
    for ii = 1:Ny
        subplot(Ny,Nu,i+(ii-1)*Nu);
        semilogx(f,20*log10(squeeze(difP(ii,i,:))));
        set(gca,'xscale','log');
        xlim([f(2) f(end)]);
        xlabel('Frequency [Hz]');
        ylabel('Estimation Error LPM [dB]');
    end
end
