function [U] = t_svd(A)
%T_SVD Summary of this function goes here
%   Detailed explanation goes here

[~,~,D]=size(A);

% if(strcmp(select,'full')==1)
%      U=zeros(H,H,D);S=zeros(H,W,D);V=zeros(W,W,D);
%     A=fft(A,[],3);
%     for i=1:D
%         [Ui,Si,Vi]=svd(A(:,:,i));
%         U(:,:,i)=Ui;S(:,:,i)=Si;V(:,:,i)=Vi;
%     end
%     
%     U=ifft(U,[],3);S=ifft(S,[],3);V=ifft(V,[],3);
%     A=ifft(A,[],3);
%     
% elseif(strcmp(select,'cut')==1)
%     A=fft(A,[],3);
%     for i=1:D
%         [Ui,Si,Vi]=svd(A(:,:,i));
%         ss=diag(Si);[~,index]=sort(ss,'descend');
%         U(:,:,i)=Ui(:,index(1:R));S(:,:,i)=Si(index(1:R),index(1:R));V(:,:,i)=Vi(:,index(1:R));
%     end
%     U=ifft(U,[],3);S=ifft(S,[],3);V=ifft(V,[],3);
%     A=ifft(A,[],3);
% end


A=fft(A,[],3);

for i=1:D
    [Ui,~,~]=svd(A(:,:,i),'econ');
    %     ss=diag(Si);[~,index]=sort(ss,'descend');
%     Si = diag(Si);
%     r = sum(Si>1e-4*Si(1));
%     Ui = Ui(:,1:r);
%     Si = Si(1:r);
%     Ui = Ui*diag(sqrt(Si));size(Ui)
%     Ui = normr(Ui);
    U(:,:,i)=Ui;
%     S(:,:,i)=Si;V(:,:,i)=Vi;
end
U=ifft(U,[],3);
% S=ifft(S,[],3);V=ifft(V,[],3);


end


