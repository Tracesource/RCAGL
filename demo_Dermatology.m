clear;
clc;
warning off;
addpath(genpath('./'));

%% dataset
ds = {'Dermatology'};
dsPath = '.\datasets\';


for dsi =1:1:length(ds)
    dataName = ds{dsi}; disp(dataName);
    load(strcat(dsPath,dataName));
    k = length(unique(Y));
    n = length(Y);
    
    lambda = [1 100 1000 10^6];
    anchor = [k,2*k,3*k];
        
    %%
    allresult = [];
    for ichor = 1:length(anchor)
        for id = 1:length(lambda)
            tic;
            [A,C,D,iter,obj,y1] = rcagl(X,Y,anchor(ichor),lambda(id));
            res = Clustering8Measure(y1, Y);
            timer(ichor,id)  = toc;
            fprintf('Anchor:%d \t Lambda:%d\t Res:%12.6f %12.6f %12.6f %12.6f \tTime:%12.6f \n',[anchor(ichor) lambda(id) res(1) res(2) res(3) res(4) timer(ichor,id)]);
            allresult = [allresult;res timer(ichor,id)];
        end
    end
    [c,d] = max(allresult(:,1));
    maxresult = allresult(d,:);
end


