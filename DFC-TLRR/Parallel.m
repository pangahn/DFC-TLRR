function [result, progress] = Parallel(dat, gnd, settings, progress)
%% initialization parameters
s = settings.s;
lambda = settings.lambda;
nClass = settings.nClass;
partition = settings.partition;

fprintf(['partition = ', num2str(partition), ', lambda = ', num2str(lambda),...
    ', s = ', num2str(s(1)), ...
    ', run = ', num2str(progress.count+1), '/', num2str(progress.total) '\n']);

%% D step
%  Divide input tensor into sub-tensor
[partitions, dat_new, gnd] = Divide(dat, gnd, partition);

%% F step
%  Factor submatrices in parallel
maxtime = 0;
[~, N, D] = size(dat);
Z_admm = zeros(N, N, D);

for i = 1:partition
    C = dat(:, partitions{i}, :);
    tic1 = tic;
    scaling_term = sqrt(N/size(C, 2));
    Z = Admm(dat_new, C, scaling_term * lambda, settings); 
    Z_admm(:, partitions{i}, :) = Z; 
    subtime = toc(tic1);
    if subtime > maxtime
        maxtime = subtime;
    end
end

%% Combine step
%  Combine submatrix estimates
tic2 = tic;
Z1 = Z_admm(:,partitions{partition},:);
Z1_f = fft(Z1,[],3);
[~,~,D] = size(Z1);
for  i = 1:D
    [Ui,~,~] = svd(Z1_f(:,:,i),'econ'); 
    Proj(:,:,i) = Ui*Ui';
end
Proj = ifft(Proj,[],3);
Z_Proj = [];
for i = 1:partition
    if(i == partition)
        Z_proj_i = Z_admm(:,partitions{partition},:);
        Z_Proj = [Z_Proj Z_proj_i];
    else
        Zi = Z_admm(:,partitions{i},:);
        Z_proj_i = tproduce(Proj,Zi);
        Z_Proj = [Z_Proj Z_proj_i];
    end
end
protime = toc(tic2);

%% reduce
% acc, nmi, time, k, lambda, nokeep
result = zeros(length(s), 6);
for nk=1:length(s)
    progress.count = progress.count + 1;
    
    if nk ~= 1
       fprintf(['partition = ', num2str(partition), ...
           ', lambda = ', num2str(lambda),...
           ', s = ', num2str(s(nk)), ...
           ', run = ', num2str(progress.count), ...
           '/', num2str(progress.total) '\n']); 
    end
  
    tic3 = tic;
    Z = reduce_Z(Z_Proj, s(nk));
    W = compute_W(Z);
    protime2 = toc(tic3);
    time = maxtime + protime + protime2;
    
    % construct W
    [U,S,~] = svd(W,'econ');
    S = diag(S);
    r = sum(S>1e-4*S(1));
    U = U(:,1:r); S = S(1:r);
    U = U*diag(sqrt(S));
    U = normr(U);
    W = (U*U').^4;
    
    % Spectral Clustering
    grps = SpectralClustering(W, nClass);
    [~, nmi, ~] = compute_nmi(gnd, grps);
    grps_best = bestMap(gnd, grps); % Very time-consuming
    acc = length(find(gnd == grps_best))/length(gnd);

%     confusion_chart(gnd, P_label); % plot confusion chart

    fprintf(['\t*** the accuracy is: ' num2str(acc) '\n']);
    fprintf(['\t*** the nmi score is: ' num2str(nmi) '\n']);
    fprintf(['\t*** the time cost is : ' num2str(time) '\n']);
    output = [acc nmi time partition lambda s(nk)];
    dlmwrite(progress.filename, output, '-append', 'precision',...
        '%.4f', 'newline', 'pc' );
    result(nk, :) = output;
end
end

function Z=reduce_Z(Z, nokeep)
if nokeep~=0
    [H, L, ~] = size(Z);
    tmp = zeros(H, 1);
    for i=1:L
        for j=1:H
            cc = Z(j, i, :);
            a = norm(cc(:));
            tmp(j) = a;
        end
        [~, idx] = sort(tmp);
        Z(idx(1:floor(nokeep*H)), i, :) = 0;
    end
end
end

function W = compute_W(Z_proj)
%COMPUTE_W Summary of this function goes here
%   Detailed explanation goes here
[H1, H2] = size(Z_proj(:,:,1));
W = zeros(H1, H2);
for j = 1:H1
    for k = 1:H2
        W(j, k) = norm(tensor(Z_proj(j,k,:)))+norm(tensor(Z_proj(k,j,:)));
    end
end

end