function [normX, method] = normalize(X, method)
%NORMALIZEFRO normalize data matrix X, 
% each column of X is a image data.
m = size(X, 1);

if strcmp(method, 'l2')
    % divide each element of X by L2-norm of its column
    normX  =  X ./ repmat(sqrt(sum(X.*X)), [m,1]);
    method = 1; % for writing
elseif strcmp(method, 'max')
    % rescale to a unit
    normX = X./255;
    method = 2;
else
    disp('unknown method!')
    normX = X;
    method = 0;
end

end