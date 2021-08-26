function [G_LPM] = LPMClosedLoopPeriodicRobustFRM(u,y,r,n,R,P)
% This script will calculate a local polynomial model for the given reference, in- and output of a system.
% The system is assumed to be SISO and in closed loop, see figure 7-4 Pintelon2012.
% Inputs:
%     u : Input to (open loop) plant in time domain
%     y : Output of plant to given input signal u
%     r : Reference to the plant (not disturbance! see Pintelon 2012 Fig. 7-4)
%     n : Window size (left and right) for frequency bin k
%     R : Degree of polynomial, e.g. G(omega_k+r) = G(omega)+sum_s=1^R g_s(k)*r^s
%     P : Amount of periods in signal r
% Outputs:
%     G_LPM : Estimated plant dynamics FRM in robust manner
%% errors
dof = 2*n+1-(R+1)*2;
if dof<1
    error(['Not high enough DOF = ',num2str(dof),' < 1']);
elseif size(u,2)>1 || size(y,2)>1
    error('Not yet defined for MIMO systems or provide in- and output data as vectors');
end
%% define variables.
Ns = length(u);  % Total amount of samples
Np = Ns/P;       % Samples in one period
N = floor(Ns/2); % amount of samples up to nyquist
Nk = floor(Np/2);

thetaHat = zeros(2,(R+1),Nk); % Pintelon2012 (7-71)
K1 = @(r) (r*ones(R+1,1)).^((0:R)'); % basis for LPM

Uf=fft(u)/sqrt(P*Np); % correct Pintelon2012 (7-66)
Yf=fft(y)/sqrt(P*Np); % correct Pintelon2012 (7-66)
Rf=fft(r)/sqrt(P*Np); % correct Pintelon2012 (7-66)

Yk = Yf(1:P:N)'; % up to nyquist frequency (!!!!???)
Uk = Uf(1:P:N)';
Zk = [Yk;Uk]; % Pintelon 2012 (7-48)
Rk = Rf(1:P:N)';



% Enzo's way
% Y(n + 1 : N + 2 * n ) = Yf(1 : N + n);
% U(n + 1 : N + 2 * n ) = Uf(1 : N + n);
% Y(1 : n )  = conj(Yf(n + 1 : -1 : 2 ));
% U(1 : n )  = conj(Uf(n + 1 : -1 : 2 ));
%% remove noise transient terms
for k = 1:Nk
% for k = n + 1 : N + n; % Enzo's way (r=-n:n by default)
    if k<n+1 % left border Pintelon2012 (7-29)
        p = n-k+1;
        r=-n+p:n+p;
    elseif k>Nk-n % right border Pintelon2012 (7-29)
        p=-n+Nk-k;
        r=-n+p:n+p;
    else % everything else
        r = -n:n;
    end
    Kn = zeros((R+1),2*n+1); % reset Kn for every iteration k
    for i = 1:2*n+1
        Kn(:,i) = K1(r(i)); % Pintelon2012 between (7-71) and (7-72)
    end
    
    % scaling, see Pintelon2012 (7-25)
    Dscale = zeros(R+1);
    for i = 1:(R+1)
        Dscale(i,i) = norm(Kn(i,:),2);
    end
    
    Kn = Dscale\Kn;
    
    [U_k,S_k,V_k] = svd(Kn'); % better computational feasability Pintelon 2012 (7-24)
    thetaHat(:,:,k) = Zk(:,k+r)*U_k/S_k'*V_k';
    thetaHat(:,:,k) = thetaHat(:,:,k)/Dscale;
end
THz_LPM = squeeze(thetaHat(:,1,:)); % calculate LPM estimate of system
Zkh=Zk-THz_LPM; % subtract noise transient from Z(k) Pintelon2012 (7-73)

G_LPM = Zkh(1,:)./Zkh(2,:);

% %% loop over frequency bins
% for k = 1:N
% % for k = n + 1 : N + n; % Enzo's way (r=-n:n by default)
%     if k<n+1 % left border Pintelon2012 (7-29)
%         p = n-k+1;
%         r=-n+p:n+p;
%     elseif k>N-n % right border Pintelon2012 (7-29)
%         p=-n+N-k;
%         r=-n+p:n+p;
%     else % everything else
%         r = -n:n;
%     end
%     Kn = zeros(2*(R+1),2*n+1); % reset Kn for every iteration k
%     for i = 1:2*n+1
%         Kn(:,i) = [K1(r(i)).*Rk(k+r(i)); K1(r(i))]; 
%     end
%     
%     % scaling, see Pintelon2012 (7-25)
%     Dscale = zeros(2*(R+1));
%     for i = 1:2*(R+1)
%         Dscale(i,i) = norm(Kn(i,:),2);
%     end
%     
%     Kn = Dscale\Kn;
%     
%     [U_k,S_k,V_k] = svd(Kn'); % better computational feasability Pintelon 2012 (7-24)
%     thetaHat(:,:,k) = Zkh(:,k+r)*U_k/S_k'*V_k';
%     thetaHat(:,:,k) = thetaHat(:,:,k)/Dscale;
% end
% Grz_LPM = squeeze(thetaHat(:,1,:)); % calculate LPM estimate of system
% Trz_LPM = squeeze(thetaHat(:,R+2,:)); % calculate LPM estimate transient contribution of system 
% G_LPM = Grz_LPM(1,:)./Grz_LPM(2,:);
% T_LPM = Trz_LPM(1,:)./Trz_LPM(2,:);
% 
% % Enzo's way
% % G = thetaHat(n+1:end,1); 
% % T = thetaHat(n+1:end,R+2);
end

