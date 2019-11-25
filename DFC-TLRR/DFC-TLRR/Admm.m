function [Z,J,S] = Admm(M, C, lambda, settings)
%ADMM This routine solves the following optimization problem,
% min |J|_*f + lambda*|S|_F^2
% s.t., C = M*Z+S
%       Z=J
% inputs:
%   M -- H*N*D data tensor
%   C -- H*l*D data sub-tensor
%   lambda -- balance parameter
%   settings:
%     method: |J|_*f (:nuclear) or |J|_1 (:F1)
%     display: display the progress (:true/:false)

% Set and parse parameters
if isfield(settings,'method')
    method = settings.method;
else
    method = 'nuclear'; % default
end
fprintf(['\t--- method is: ' method '\n']);

if isfield(settings,'display')
    display = settings.display;
else
    display = false;
end

%% Initializing variables
[H,N,D] = size(M);
L = size(C,2);
mumax = 1e10;
maxIter = 500;
costTrack = zeros(maxIter, 1); % track the cost
errorTrack = zeros(maxIter, 1); %#ok<*NASGU> % track the error

J = zeros(N,L,D);
Z = zeros(N,L,D);
S = zeros(H,L,D);

G1 = zeros(H,L,D);
G2 = zeros(N,L,D);

mu = 0.1;
p = 1.9;
e = 1e-5;

MZ = zeros(H,L,D);
m = fft(M,[],3);
I = eye(N);

%% Start main loop

iter = 0;
while iter<maxIter
    iter = iter + 1;

    %update J
    B = Z+G2/mu; % G2 = zeros(N,L,D);
    for j = 1:D
        temp = B(:,:,j);
        if strcmp(method,'nuclear') %take F-nuclear norm
            [U,E,V] = svd(temp,'econ');
            E = sign(E).*max(abs(E)-1/mu,0);
            J(:,:,j) = U*E*V';
        elseif strcmp(method,'L1') %take l1 norm
            J(:,:,j) = sign(temp).*max(abs(temp)-1/mu,0);
        else
            error('method does not support!');
        end
    end
%     clear B temp U E V;
        
    %udpate Z
    R = C-S+G1/mu;
    R_f = fft(R,[],3); % G1 = zeros(H,L,D);
    T = fft(J-G2/mu,[],3); % G2 = zeros(N,L,D);
    z = zeros(N,L,D);
    for j = 1:D
        z(:,:,j) = (m(:,:,j)'*m(:,:,j)+I) \ (m(:,:,j)'*R_f(:,:,j) + T(:,:,j));
    end
    Z = ifft(z,[],3);
%     clear R R_f T z
    
    %update S
    z = fft(Z,[],3); % Z = zeros(N,L,D);
    for j = 1:D % M -- H*N*D
        MZ(:,:,j) = m(:,:,j) * z(:,:,j);
    end
    MZ = ifft(MZ,[],3); % MZ = zeros(H,L,D);
    S = mu / (2 * lambda + mu) .* (C - MZ + G1 ./ mu);
%     clear z;
    
    % track the cost value
    % costTrack(iter) = evalCost(J, S, lambda); 

    stop = max(max(max(max(abs(C-MZ-S)))),max(max(max(abs(Z-J)))));
    % errorTrack(iter) = stop;
    if display && (iter==1 || mod(iter,50)==0 || stop<e)
        fprintf(['\t+++ iter = ' num2str(iter) ...
            ', mu=' num2str(mu,'%2.1e')...
            ', stopALM=' num2str(stop,'%2.3e') '\n']);
    end
    
    if stop<e
        % plotTrack(costTrack, 'cost value')
        % plotTrack(errorTrack, 'error')
        break;
    else
        %update G1,G2,mu
        G1 = G1 + mu * (C - MZ - S);
        G2 = G2 + mu * (Z - J);
        mu = min(p * mu, mumax);
    end
end

end

function cost = evalCost(J, S, lambda)
    [~,~,l]  = size(J);
    nulear_norm = 0;
    for i=1:l
        nulear_norm = nulear_norm + norm_nuclear(J(:,:,i));
    end
    cost = nulear_norm + lambda * sqrt(sum(S(:).^2));
end

function plotTrack(track, string)
    lastNoneZero = find(track~=0, 1, 'last');
    figure
    plot(track(1:lastNoneZero), 'r-', 'linewidth',1.5)
    xlabel('iteration')
    ylabel(string)
end
