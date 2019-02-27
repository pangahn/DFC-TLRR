function [partitions, dat, gnd] = Divide(dat, gnd, t)
%DIVIDE split a tensor dat into roughly equal-sized pieces.
%   INPUT
%     dat: A tensor dat, each lateral slice is a image
%     gnd: label of image
%     t: number of segmentation(s)
%   OUTPUT
%     partitions: cells type date, each cell is a index list
%     dat: data after shuffling
%     gnd: labels after shuffling

N = size(dat,2);
seq = randperm(N);
dat = dat(:, seq, :);
gnd = reshape(gnd(seq),[],1);
partitions = cell(t, 1);
splitSize = 1/t * N;

for i=1:t
   idx = round((i-1)*splitSize + 1):round(i*splitSize);
   
%    %下面代码用于统计每类样本在每个分割中的分布数据 
%    stat = tabulate(gnd(idx));
%    dlmwrite('Stat.csv', stat, '-append', 'precision', '%.4f', 'newline', 'pc' );
%    fp = fopen('Stat.csv','a');
%    fprintf(fp, '\n');
%    fclose(fp);

   partitions{i} = seq(idx);
end
end