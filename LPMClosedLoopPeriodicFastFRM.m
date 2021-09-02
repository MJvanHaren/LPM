function [G_LPM] = LPMClosedLoopPeriodicFastFRM(u,y,r,n,R,P)
% This script will calculate a local polynomial model for the given reference, in- and output of a system.
% The system is assumed to be in closed loop, see figure 7-4 Pintelon2012.
% Inputs:
%     u : Input to (open loop) plant in time domain
%     y : Output of plant to given input signal u
%     r : Reference to the plant (not disturbance! see Pintelon 2012 Fig. 7-4)
%     n : Window size (left and right) for frequency bin k
%     R : Degree of polynomial, e.g. G(omega_k+r) = G(omega)+sum_s=1^R g_s(k)*r^s
% Outputs:
%     G_LPM : Estimated plant dynamics
%     T_LPM : Estimated transient dynamics
%% errors
dof = 2*n+1-(R+1)*2;
if dof<1
    error(['Not high enough DOF = ',num2str(dof),' < 1']);
end
%% define variables.
N = length(u);      % Total amount of samples
Np = N/P;           % Samples in one period
Nn = floor(Np/2);   % amount of samples per period up to nyquist

Nu = size(u,2); % number of inputs
Ny = size(y,2); % number of outputs

thetaHat = zeros(2*Ny,(Nu+1)*(R+1),Nn); % Pintelon2012 (7-6)
K1 = @(r) (r*ones(R+1,1)).^((0:R)'); % basis for LPM

Uf=fft(u)/sqrt(P*Np); % correct Pintelon2012 (7-66)
Yf=fft(y)/sqrt(P*Np); % correct Pintelon2012 (7-66)
Rf=fft(r)/sqrt(P*Np); % correct Pintelon2012 (7-66)

Yk = Yf(1:P:P*Nn,:)'; % up to nyquist frequency (!!!!???) ESSENTIAL
Uk = Uf(1:P:P*Nn,:)';
Zk = [Yk;Uk];       % Pintelon 2012 (7-48)
Rk = Rf(1:P:P*Nn,:)';

%% loop over frequency bins
Grz_LPM = zeros(2*Ny,Nu,Nn);
G_LPM = zeros(Ny,Nu,Nn);
for k = 1:Nn
    if k<n+1 % left border Pintelon2012 (7-29)
        p = n-k+1;
        r=-n+p:n+p;
    elseif k>Nn-n % right border Pintelon2012 (7-29)
        p=-n+Nn-k;
        r=-n+p:n+p;
    else % everything else
        r = -n:n;
    end
    Kn = zeros((1+Nu)*(R+1),2*n+1); % reset Kn for every iteration k
    for i = 1:2*n+1
        Kn(:,i) = [kron(K1(r(i)),Rk(:,k+r(i))); K1(r(i))]; 
    end
    
    % scaling, see Pintelon2012 (7-25)
    Dscale = zeros((1+Nu)*(R+1));
    for i = 1:(1+Nu)*(R+1)
        Dscale(i,i) = norm(Kn(i,:),2);
    end
    
    Kn = Dscale\Kn;
    
    [U_k,S_k,V_k] = svd(Kn'); % better computational feasability Pintelon 2012 (7-24)
    thetaHat(:,:,k) = Zk(:,k+r)*U_k/S_k'*V_k';
    thetaHat(:,:,k) = thetaHat(:,:,k)/Dscale;
    Grz_LPM(:,:,k) = thetaHat(:,:,k)*[eye(Nu);zeros(size(thetaHat(:,:,k),2)-Nu,Nu)];% calculate LPM estimate of system
    G_LPM(:,:,k) = Grz_LPM(1:Nu,:,k)./Grz_LPM(Nu+1:end,:,k);
end
end

