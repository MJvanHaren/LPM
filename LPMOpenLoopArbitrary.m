function [G_LPM, T_LPM] = LPMOpenLoopArbitrary(u,y,n,R)

%% errors
dof = 2*n+1-(R+1)*2;
if dof<1
    error(['Not high enough DOF = ',num2str(dof),' < 1']);
end
if size(u,2)>1 || size(y,2)>1
    error('Not yet defined for MIMO systems')
end
%% define variables.
Np = length(u);
N = floor(Np/2);
thetaHat = zeros(N,2*(R+1));
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
G_LPM = thetaHat(:,1);
T_LPM = thetaHat(:,R+2);

% Enzo's way
% G = thetaHat(n+1:end,1); 
% T = thetaHat(n+1:end,R+2);
end

