function main(parms)
%MAIN is the run script of experiment

%% initialization parameters
dataSet = parms.dataSet;
lambda = parms.lambda;
normalizeMethod = parms.normalizeMethod;
partition = parms.partition;
run = parms.run;
s = parms.s;

%% data processing
fprintf(['===== ' dataSet ' =====\n'])
load (dataSet);   %#ok<LOAD> <dat, gnd>

if contains(dataSet,'MNIST_Test')
    dat = T_test_images;
    gnd = test_labels;
    clear T_test_images test_labels;
elseif contains(dataSet,'MNIST_Train')
    dat = T_train_images;
    gnd = train_labels;
    clear T_train_images train_labels;
end

nClass = length(unique(gnd));
dat = Normalize(dat, normalizeMethod); %#ok<*NODEF>

%% setting up result folds
filename = fullfile('Output', [dataSet ,'.csv']);
if ~exist(filename, 'file')
    header = {'accuracy', 'nmi', 'time', 'partition', 'lambda', 's', 'loop'};
    fid = fopen(filename,'a+');
    fprintf(fid,'%s,%s,%s,%s,%s,%s,%s\n',header{1},header{2},header{3},...
        header{4},header{5},header{6},header{7});
    fclose(fid);
end

filename2 = fullfile('Output', [dataSet ,'_mean.csv']);
if ~exist(filename2,'file')
    header = {'mean_acc', 'mean_nmi', 'std_acc', 'std_nmi', 'partition', ...
        'lambda', 's', 'time'};
    fid = fopen(filename2,'a+');
    fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s\n',header{1},header{2},...
        header{3},header{4},header{5}, header{6}, header{7}, header{8});
    fclose(fid);
end

%% main loops
progress.total = length(lambda)*length(partition)*length(s)*run;
progress.count = 0;
progress.filename = filename;

settings.s = s;
settings.nClass = nClass;
settings.method = parms.method;
settings.display = parms.display;

for j=1:length(partition)
    for i=1:length(lambda)
        % accuracy, nmi, time, partition, lambda, s
        summary = zeros(length(s), 6, run);
        
        settings.lambda = lambda(i);
        settings.partition = partition(j);
        for iter = 1:run
            progress.iter = iter;
            [result, progress] = Parallel(dat, gnd, settings, progress);
            summary(:, :, iter) = result;
        end
        for kk=1:length(s)
            acc_mean = mean(summary(kk, 1, :));
            nmi_mean = mean(summary(kk, 2, :));
            std_acc = std(summary(kk, 1, :));
            std_nmi = std(summary(kk, 2, :));
            time_mean = mean(summary(kk, 3, :));
            output = [acc_mean nmi_mean std_acc std_nmi partition(j) ...
                lambda(i) summary(kk, 6, 1) time_mean];
            dlmwrite(filename2, output, '-append', 'precision',...
                '%.4f', 'newline', 'pc' );
        end
    end
end
end