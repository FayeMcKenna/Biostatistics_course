%class 5
%clustering

%% clean up
clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same
% one dimension shows the most variation in one way and the other dimension
% the least variation, look at broad patterns and use those
%load pancreas data set
try
    A;
catch
    load melanoma_working_data.mat
end    


% correlation matrix , 
% hierarchical clustering is finding the 'closeness' as a way to order in
% order to put closer things together for the matrix
% combine the two things that are closest together then recombine them to
% create a new matrix tree
% continue to find closest and average

A_log=log10(1+A);
figure; imagesc(A_log);

%how many genes are not expressed in any cell?

length(find(sum(A')==0))

% condition_1=expressed (give top half)
condition_1=find(sum(A')>median(sum(A'))); length(condition_1)
% condition_2=standard deviation 
condition_2=find(std(A')>2);

% condition_3= number of cells expressed in , take the matrix whatever is
% higher than one sum in in one way
B= A_log; B(find(B>0))=1;


condition_3= find(sum(B')>500);

dyn_genes= intersect(condition_3,intersect(condition_1,condition_2));length(dyn_genes)

figure;imagesc(A_log(dyn_genes,:))


C=1-squareform(pdist(A_log(dyn_genes,:)', 'correlation')); % added correlation and 1- to find correaltion not euclidian

figure;imagesc(C);

% now hclustering

L = linkage(C)

% [H,T,OUTPERM= dendrogram(....)

[H,T,OUTPERM]= dendrogram(L,0); % for every cell what order it should be only be outperm
    

linkage(C); % creates a hcluster
figure;imagesc(C(OUTPERM,OUTPERM));

% the diaganol line should be blue distance from itself is zero because it
% is euclidian distance

% maybe we want to use correlation?

% check what is each gene in matrix?

% k mean clustering, means you have to give it the tree and the parameter
% K, how many clusters you have 
% randomly puts down 'k' centroids and the data clusters around point, then
% run again euristic algorithm(always different answer) eventually converges, Hclustering is
% deterministic
km_param=50;
km=kmeans(A_log(dyn_genes,:)',km_param);

% change order according to k

[i,xi]=sort(km);
% figure;imagesc(C)

C= corr(A_log(dyn_genes,:));
figure;imagesc(C(xi,xi)); %figure;imagesc(C(xi,xi), [0.2 0.5]);


























