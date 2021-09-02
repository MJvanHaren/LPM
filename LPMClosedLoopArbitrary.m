function [G_LPM] = LPMClosedLoopArbitrary(u,y,r,n,R)
% This script will calculate a local polynomial model for the given reference, in- and output of a system.
% The system is assumed to be SISO and in closed loop, see figure 7-4 Pintelon2012.
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
Np = length(u); % signal length
N = floor(Np/2); % amount of samples up to nyquist

Nu = size(u,2); % number of inputs
Ny = size(y,2); % number of outputs

%% Transient example from Voorhoeve2018
s=tf('s');
T = 100/(s^2+2*0.1*(15*2*pi)*s+(15*2*pi)^2);                    % Voorhoeve2018 Figure 2.2
f = linspace(0, 1 - 1/N, N) * (1/250)/2;                        % available frequencies (Fix Ts)
Tk = squeeze(frd(c2d(T,1/250,'tustin'),f).ResponseData);        % (Fix Ts)
%% dft etc
thetaHat = zeros(2*Ny,(Nu+1)*(R+1),N); % Pintelon2012 (7-6)
K1 = @(r) (r*ones(R+1,1)).^((0:R)'); % basis for LPM

Uf=fft(u)/sqrt(Np);
Yf=fft(y)/sqrt(Np);
Rf=fft(r)/sqrt(Np);
Yk = Yf(1:N,:)'; % up to nyquist frequency
Uk = Uf(1:N,:)';
Zk = [Yk;Uk]; % Pintelon 2012 (7-48)
Rk = Rf(1:N,:)';
%% loop over frequency bins
Grz_LPM = zeros(2*Ny,Nu,N);
G_LPM = zeros(Ny,Nu,N);
for k = 1:N
    if k<n+1 % left border Pintelon2012 (7-29)
        p = n-k+1;
        r=-n+p:n+p;
    elseif k>N-n % right border Pintelon2012 (7-29)
        p=-n+N-k;
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
    % T_LPM = zoveelste rij van thetahat
end


end

