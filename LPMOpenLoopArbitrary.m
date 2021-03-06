function [G_LPM, T_LPM] = LPMOpenLoopArbitrary(u,y,n,R)
% This script will calculate a local polynomial model for the given in- and output of a system.
% The system is assumed to be SISO and in open loop.
% Inputs:
%     u : Input to (open loop) plant in time domain
%     y : Output of plant to given input signal u
%     n : Window size (left and right) for frequency bin k
%     R : Degree of polynomial, e.g. G(omega_k+r) = G(omega)+sum_s=1^R g_s(k)*r^s
% Outputs:
%     G_LPM : Estimated plant dynamics
%     T_LPM : Estimated transient dynamics
%% errors
dof = 2*n+1-(R+1)*2;
if dof<1
    error(['Not high enough DOF = ',num2str(dof),' < 1']);
elseif size(u,2)>1 || size(y,2)>1
    error('Not yet defined for MIMO systems or provide in- and output data as vectors');
end
%% define variables.
Np = length(u); % signal length
N = floor(Np/2); % amount of samples up to nyquist

thetaHat = zeros(N,2*(R+1)); % Pintelon2012 (7-6)
K1 = @(r) (r*ones(R+1,1)).^((0:R)'); % basis for LPM

Uf=fft(u)/sqrt(Np);
Yf=fft(y)/sqrt(Np);
Y = Yf(1:N)'; % up to nyquist frequency
U = Uf(1:N)';

% Enzo's way
% Y(n + 1 : N + 2 * n ) = Yf(1 : N + n);
% U(n + 1 : N + 2 * n ) = Uf(1 : N + n);
% Y(1 : n )  = conj(Yf(n + 1 : -1 : 2 ));
% U(1 : n )  = conj(Uf(n + 1 : -1 : 2 ));

%% loop over frequency bins
for k = 1:N
% for k = n + 1 : N + n; % Enzo's way (r=-n:n by default)
    if k<n+1 % left border Pintelon2012 (7-29)
        p = n-k+1;
        r=-n+p:n+p;
    elseif k>N-n % right border Pintelon2012 (7-29)
        p=-n+N-k;
        r=-n+p:n+p;
    else % everything else
        r = -n:n;
    end
    Kn = zeros(2*(R+1),2*n+1); % reset Kn for every iteration k
    for i = 1:2*n+1
        Kn(:,i) = [K1(r(i)).*U(k+r(i)); K1(r(i))]; 
    end
    
    % scaling, see Pintelon2012 (7-25)
    Dscale = zeros(2*(R+1));
    for i = 1:2*(R+1)
        Dscale(i,i) = norm(Kn(i,:),2);
    end
    
    Kn = Dscale\Kn;
    
    [Uk,Sk,Vk] = svd(Kn'); % better computational feasability Pintelon 2012 (7-24)
    thetaHat(k,:) = Y(k+r)*Uk/Sk'*Vk';
    thetaHat(k,:) = thetaHat(k,:)/Dscale;
    Vhn = Y(k+r)'-thetaHat(k,:)*Kn; % Not used right now, Pintelon (7-11)
end
G_LPM = thetaHat(:,1); % calculate LPM estimate of system
T_LPM = thetaHat(:,R+2); % calculate LPM estimate transient contribution of system 

% Enzo's way
% G = thetaHat(n+1:end,1); 
% T = thetaHat(n+1:end,R+2);
end

