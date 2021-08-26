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
G0s = C/(1+P*C);
G0 = c2d(G0s,Ts,'tustin');
ssG0 = ss(G0);

n = 6 ; % window size left and right
K1 = @(r) [r.^0;r.^1;r.^2;r.^3;r.^4]; % basis for LPM

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


%% simulate
y = lsim(G0,u,t);

Uf=fft(u)/sqrt(Np);
Yf=fft(y)/sqrt(Np);

% Select the relevant frequency range
N = floor(Np/2);
Y = Yf(1:N)';
U = Uf(1:N)';

% Enzo's way
% Y(n + 1 : N + 2 * n ) = Yf(1 : N + n);
% U(n + 1 : N + 2 * n ) = Uf(1 : N + n);
% Y(1 : n )  = conj(Yf(n + 1 : -1 : 2 ));
% U(1 : n )  = conj(Uf(n + 1 : -1 : 2 ));


Getfe = etfe([y u],500,N);
%% LPM
R = length(K1(1))-1;
dof = 2*n+1-(R+1)*2;
if dof<1
    error(['DOF not high enough = ',num2str(dof)]);
end
thetaHat = zeros(N,2*(R+1));

for k = 1:N
% for k = n + 1 : N + n; % Enzo's way (r=-n:n by default)
    if k<n+1 % left border
        p = n-k+1;
        r=-n+p:n+p;
    elseif k>N-n % right border
        p=-n+N-k;
        r=-n+p:n+p;
    else % everything else
        r = -n:n;
    end
    Kn = zeros(2*(R+1),2*n+1);
    for i = 1:2*n+1
        Kn(:,i) = [K1(r(i)).*U(k+r(i)); K1(r(i))]; 
    end
    
    % scaling, see Pintelon2012 (7-25)
    Dscale = zeros(2*(R+1));
    for i = 1:2*(R+1)
        Dscale(i,i) = norm(Kn(i,:),2);
    end
    
    Kn = Dscale\Kn;
    
    [Uk,Sk,Vk] = svd(Kn');
    thetaHat(k,:) = Y(k+r)*Uk/Sk'*Vk';
    thetaHat(k,:) = thetaHat(k,:)/Dscale;
    Vhn = Y(k+r)'-thetaHat(k,:)*Kn;
end
G = thetaHat(:,1);
T = thetaHat(:,R+2);
% G = thetaHat(n+1:end,1); % Enzo's way
% T = thetaHat(n+1:end,R+2);
%%
figure
opts = bodeoptions;
opts.FreqUnits = 'Hz';
opts.xlim = [1e-1 125];
bodemag(G0,opts); hold on;
semilogx(f,20*log10(abs(G)),'Color',c2); hold on %-(n-1)*Ts or -n*Ts is perfect! -> fix in freq window?
semilogx(f,20*log10(abs(T(1:length(f)))),'Color',c3);
semilogx(Getfe.Frequency*128/pi,20*log10(squeeze(abs(Getfe.ResponseData))),'Color',c4);
legend('True plant','Estimated plant','Estimated transient plant','ETFE')

figure
[magG,~,~] = bode(G0,f*2*pi);
semilogx(f,20*log10(abs(G))-squeeze(magG),'Color',c2); hold on 
[magG,~,~] = bode(G0,Getfe.Frequency*128*2);
semilogx(Getfe.Frequency*128/pi,20*log10(squeeze(abs(Getfe.ResponseData)))-squeeze(magG),'Color',c3);
xlabel('Frequency [Hz]'); xlim([f(1) f(end)]);
ylabel('Estimation Error [dB]');



