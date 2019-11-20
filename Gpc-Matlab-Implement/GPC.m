function GPC
clear all;clf;clc;
N      = 3;                   % Preview Horizon
Nu     = 3;                   % Control Horizon
A      = [1, -0.9, 0.6, 0.3];
B      = [1, 0.5];            % Order of B
na     = length(A) - 1;       % Order of A
nb     = length(B) - 1;
noise  = 0.1;                 % Amplitude of noise
alpha  = 0.2;                 % Smooth Factor
rho    = 0.99;                % Forgotten Factor
lambda = 0.1;                 % Weight Gain for intputs
nIterationCount = 500;
% Yr              = zeros(1,nIterationCount+100);
% Yr(20:end)      = 10;
Yr              = sin((1:nIterationCount+100)/30);
initOutput      = zeros(size(A));
initInput       = zeros(size(B));

% Initial Covariance
P=eye(na+nb+1);		% covariance of parameters

% hatU    = [0 initInput];
hatU    = zeros(nb + 2, 1);
Y       = zeros(na + 1, 1);
theta   = ones(na+nb+1, 1);
% theta = [hatA; hatB];


out   = zeros(1,nIterationCount);
in    = zeros(1,nIterationCount);

for count = 1 : nIterationCount
	% Step 1: Calculate the actual output according to the general model
	Phi = flipud([diff(hatU); -diff(Y)]);
	y   = ((-1)*A(2:end))*Y(na:-1:1)+B*hatU(end:-1:2)+(rand(1)-0.5)*2*noise;
	Y   = [Y(2:end); y];
	out(count) = Y(na + 1);
	deltaY     = Y(end)-Y(end-1);


	theta = theta + (rho + Phi' * P * Phi)\ P * Phi * (deltaY - Phi' * theta);
	P     = (P-(rho+Phi'*P*Phi)\P*Phi*Phi'*P)/rho;

	hatA    = [1 theta(1:na)'];		% estimation of coeffeciences of A
	hatB    = theta(na+1:end)';		% estimation of coeffeciences of B

	% ------------------------------------
	% Estimate E,F,G,H using diophantine iteration
	% ------------------------------------
	% Solve the Diophantine eqaution
	% $1      = E_{j}A\Delta+BF_{j}$
	% $E_{j}B = G_{j}+z^{-j}H_{J}$
    [F, H, g] = Diophantine(conv([1,-1],hatA), hatB, N);
	Y0 = F * Y(end:-1:1)+H* diff(hatU(2:end));

	GG=zeros(N,Nu);
	for n=1:Nu
		GG(n:N,n)=g(1:N-n+1);
	end

	Gain      = (GG'*GG+lambda*eye(Nu))\GG';		% actually, we need only the first line of the Gain matrix, thus
	Gain      = Gain(1,:);

	% Generate the smoothed trajectory
	temYr = Yr(count);
	Yd = (Y(size(hatA,2)) - temYr) * cumprod( alpha * ones(N, 1)) + temYr;
	% Calculate the control derivative
	deltaU  = Gain*(Yd-Y0);
	% Update all the matrix
    hatU       = [hatU(2:end); hatU(end) + deltaU];
	in(count)  = hatU(end);
end

figure(1);clf;
subplot(211); plot(out);
subplot(212); plot(in);

end

function [F, H, g] = Diophantine(A, B, N)
    % Since in Simulink, dynamic expand the matrix is not supported
    % Preallocate space
    e = zeros(1, N);
    g = zeros(1, N);
    F = zeros(N, length(A) - 1);
    H = zeros(N, length(B) - 1);
    % assign values
    e(1)  = 1;
	f  = -A(2:end);
 	g(1)  = e(1)*B(1);
 	h  = B(2:end);
    
	H(1,:) = h;
	F(1,:) = f;

	% Prediction for n steps
	for n=2:N
        ej = f(1);
        % $E_{j+1} = E_j + e_j z^{-j}$
        e(n) = ej;
		% $F_{j+1}=z[F_{j}-e_{j}A]$
		f = [f(2:end), 0]  - A(2:end) * ej;
        g(n) = ej*B(1)+h(1);
		h = [h(2:end), 0] + ej * B(2:end);
        
		F(n,:) = f;
		H(n,:) = h;
    end
end
