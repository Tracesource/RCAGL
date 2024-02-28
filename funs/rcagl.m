function [A,C,D,iter,obj,y1] = rcagl(X,Y,numanchor,lambda)
% m      : the number of anchor. the size of Z is m*n.
% X      : n*di

%% initialize
maxIter = 50 ; % the number of iterations
IterMax = 50;

m = numanchor;
numclass = length(unique(Y));
numview = length(X);
numsample = size(Y,1);

C = zeros(m,numsample); % m  * n
XX = [];
for p = 1 : numview
    X{p} = mapstd(X{p}',0,1);
    XX = [XX;X{p}];
end
[XU,~,~]=svds(XX',m);
rand('twister',12);
[IDX,~] = kmeans(XU,m, 'MaxIter',100,'Replicates',10);
for i = 1:numsample
    C(IDX(i),i) = 1;
end

for i = 1:numview
   di = size(X{i},2); 
   D{i} = C;
   A{i} = zeros(di,m);
end

flag = 1;
iter = 0;
%%
while flag
    iter = iter + 1;
    
    %% optimize Ai
    parfor iv=1:numview
        G = X{iv}*(C+D{iv})';      
        [U,~,V] = svd(G,'econ');
        A{iv} = U*V';
    end
    
    %% optimize Di
    parfor iv=1:numview
        H = 2*(X{iv}-0.5*A{iv}*C)'*A{iv};      
        for ii=1:numsample
            ut = H(ii,:)./(1+4*lambda);
            D{iv}(:,ii) = EProjSimplex_new(ut');
        end
    end
    
    %% optimize C
    B = zeros(numsample,m);
    G = 0;
    for iv=1:numview
        B = B+(X{iv}-0.5*A{iv}*D{iv})'*A{iv};
        G = G+D{iv}./numview;
    end
    B = B./numview;
    [P,~,~,y1] = coclustering_bipartite_fast_re(B, G', numclass,IterMax);
    C = P';

    term1 = 0;
    term2 = 0;
    for iv = 1:numview
        term1 = term1 + norm(X{iv}-0.5*A{iv}*(C+D{iv}),'fro')^2;
        term2 = term2 + norm(D{iv},'fro')^2;
    end
    
    obj(iter) = term1+lambda*term2;
    
	if (iter>1) && (abs((obj(iter-1)-obj(iter))/(obj(iter-1)))<1e-4 || iter>maxIter || obj(iter) < 1e-10)
        flag = 0;
    end
end
         
         
    
