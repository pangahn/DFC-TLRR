function[Z,E,time,gnd] = run_lrr_dfc(X, lambda, part, gnd)

%maxtime record the time of the longest running subproblem
maxtime = 0;
%protime record the time required to combine submatrix estimates 
protime = 0;

% partition columns
[partitions,gnd, X_new] = partition(1:size(X,2), part, gnd,X); 

% run lrr on each subproblem
Z = [];

for i = 1:part
  C = X(:, partitions{i});
  
  tic1 = tic;
  scaling_term = sqrt(size(X_new,2)/size(C,2));
  % rescale lambda??
  [C_hat, ~] = run_lrr(C, X_new, scaling_term * lambda);
  subtime = toc(tic1);
  if subtime > maxtime
      maxtime = subtime;
  end
  tic2 = tic;

  if i == 1
    [U,~,~] = svd(C_hat,'econ');
    proj = U * U';
  else
    C_hat = proj * C_hat;
  end
  Z = [Z C_hat];
  
  protime = protime+toc(tic2);
end
E = X_new - X_new*Z;

time = maxtime+protime;

end

function[Z,E] = run_lrr(X, A, lambda)
%%% for efficiency 
% [~,S,Q] = svd(A,'econ');
% S = diag(S);
% r = sum(S>1e-4*S(1));
% Q = Q(:,1:r);
% B = A*Q;
% [Z,E] = lrra(X,B,lambda,false);
% Z = Q * Z;

Q = orth(A');
Dict = A*Q;
[Z,E] = lrra(X, Dict, lambda, false);
Z = Q * Z;
end

function [res,gnd,X_new] = partition(v,k,gnd,X)
% PARTITION Partitions a vector v into k contiguous pieces of nearly
% equal length.

N = size(v,2);
v = randperm(N);
gnd = gnd(v);
X_new = X(:,v);

num_per_part = floor(length(v)/k);
num_long_parts = length(v) - k*num_per_part;
res = {};
part_start = 1;
for ii = 1:k
   if ii <= num_long_parts
      part_len = num_per_part+1;
   else
      part_len = num_per_part;
   end
   res{ii} = v(part_start:(part_start+part_len-1));
   part_start = part_start+part_len;
end
end
