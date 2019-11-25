function data = Normalize(data, method)
%NORMALIZEFRO normalize 3d tensor data, 
% each lateral slice of data is a image.
%
% Copyright @ Gan Phua, 2018

if strcmp(method, 'max')
    % scale the features to [0, 1]
    data = data./255;
elseif strcmp(method, 'l2')
    [~, nsample, ~] = size(data);
    for i = 1:nsample
        a = data(:,i,:);
        b = tenmat(a,1);
        data(:,i,:) = data(:,i,:) ./ max(1e-12,norm(b.data,2));
    end
end

end