function main_dfc_lrr(par)

%% initialization parameters
dataSet = par.dataSet;
lambda = par.lambda;
part = par.part;
run = par.run;
post = par.post;
normalizeMethod = par.normalizeMethod;

%% data processing
fprintf(['===== ' dataSet ' =====\n'])
load(dataSet); %#ok<*LOAD>
if contains(dataSet,'MNIST_Test')
    dat = V_test_images;
    gnd = test_labels;
    clear V_test_images test_labels;
elseif contains(dataSet,'MNIST_Train')
    dat = V_train_images;
    gnd = train_labels;
    clear V_train_images train_labels;
end
[dat, ~] = normalize(dat, normalizeMethod);
progress = length(lambda)*length(part)*run;
nClass = length(unique(gnd));
count = 0;

%% setting up result folds
filename = fullfile('Output', [dataSet ,'.csv']);
if ~exist(filename,'file')
    header = {'accuracy', 'nmi', 'time', 'partition', 'lambda', 'loop'};
    fid = fopen(filename,'a+');
    fprintf(fid,'%s,%s,%s,%s,%s,%s\n',header{1},header{2},...
        header{3},header{4},header{5},header{6});
    fclose(fid);
end

filename2 = fullfile('Output', [dataSet ,'_mean.csv']);
if ~exist(filename2,'file')
    header = {'mean_acc', 'mean_nmi', 'std_acc', 'std_nmi', ...
        'partition', 'lambda', 'time'};
    fid = fopen(filename2,'a+');
    fprintf(fid,'%s,%s,%s,%s,%s,%s,%s\n',header{1},header{2},...
        header{3},header{4},header{5},header{6},header{7});
    fclose(fid);
end

%% main loops
for j =1:length(part) 
    for i = 1:length(lambda)        
        for iter = 1:run
            % run DFC_LRR
            count = count + 1;
            
            disp(['part = ',num2str(part(j)),...
                ', lambda = ', num2str(lambda(i)),...
                ', run=', num2str(count), '/', num2str(progress)]);
            
            [Z,~,Time,labels_new] = run_lrr_dfc(dat, lambda(i), part(j), gnd);
            tic;
            tstart = tic;
            % post processing
            if post == 1
                [U,S,~] = svd(Z,'econ');
                S = diag(S);
                r = sum(S>1e-4*S(1));
                U = U(:,1:r); S = S(1:r);
                U = U*diag(sqrt(S));
                U = normr(U);
                CKSym = (U*U').^4;
            else
                CKSym = BuildAdjacency(Z, 0);
            end
            
            grps = SpectralClustering(CKSym,nClass);
            clear Z CKSym;
            
            time(j,i,iter) = Time + toc(tstart);
            
            P_label = bestMap(labels_new,grps);
            gndd = reshape(labels_new,[],1);

            accuracy(j,i,iter) = length(find(gndd == P_label))/length(gndd);
            [~, nmi(j,i,iter), ~] = compute_nmi(gndd, P_label);
            
            fprintf(['\t*** the accuracy is: ' num2str(accuracy(j,i,iter)) '\n']);
            fprintf(['\t*** the NMI score is: ' num2str(nmi(j,i,iter)) '\n']);
            fprintf(['\t*** the time cost is : ' num2str(time(j,i,iter)) '\n']);
                       
            output = [accuracy(j,i,iter) nmi(j,i,iter) time(j,i,iter)...
                part(j) lambda(i) iter];
            dlmwrite(filename, output, '-append', 'precision',...
                '%.4f', 'newline', 'pc' );
        end   
  
        std_acc(j,i) = std(accuracy(j,i,:));
        std_nmi(j,i) = std(nmi(j,i,:));
        acc_mean(j,i) = mean(accuracy(j,i,:));
        nmi_mean(j,i) = mean(nmi(j,i,:));        
        time_mean(j,i) = mean(time(j,i,:));
        
        output = [acc_mean(j,i) nmi_mean(j,i) std_acc(j,i) ...
            std_nmi(j,i) part(j) lambda(i) time_mean(j,i)];
        dlmwrite(filename2, output, '-append', 'precision',...
            '%.4f', 'newline', 'pc' );
    end    
end
end