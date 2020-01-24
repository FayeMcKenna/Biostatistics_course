
%% clean up %Hmk3
clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same

load C_elegans_dataset.mat
load('CE_genes.mat')

A= C_elegans_RNA_Seq_Embryo;
%Calculate number of genes for each cell
A_log= log10(1+A);
index=A_log~=0;
Num_of_genes=sum(index,1);

%Question A)

% Identify 500 dynamic genes in the C_elegans dataset
set1= find(std(A')> 2); %want standard deviation larger than 2
set2= find(sum(A')>3000); %sum of expression of cells >3000
set1_2= intersect(set1, set2);
B= log10(1+A); %log scaled vector of original A values
dynamic_genes= B(set1_2,:); % gives 677 genes
% now let's sort by total expression and choose top 500
[s, i]=sort(sum(dynamic_genes'), 'descend');
dyn_genes= dynamic_genes(i(1:500),:);

[~, xi ] = sort(dyn_genes);
%display the gene names of the genes with the highest coefficients
gene_names=CE_genes(xi(1:500))

figure;
imagesc(dyn_genes);
colorbar;
xlabel('Cell Number');
ylabel('Dynamic Gene Number');
title('Final Dynamic Genes Unsorted Heatmap');

% Calculate pairwise correlations of these dynamic genes
B= log10(1+dyn_genes);
% gene_i = B(500,CE_genes);
x = B();
y = B();
C = corr(x',y');
C(find(isnan(C))) = 0;
[~, xi] = sort(C,'descend');
D = xi(1:500)
C(xi(1:10)) % find correlation       

imagesc(C);

figure; corrplot(B)

%identify the most positively correlated pair and most negatively correlated pair in this gene set. 
maxC = max(max(triu(C,1))) % highest pos correlation
[row,col] = find(C==maxC,1,'first')

%     0.9857 row = 428 , col = 406 :

gene_i = gene_names(428,1) %'csn-2'
gene_j = gene_names(406,1) %'sra-35'


minC = min(min(triu(C,1))) % highest neg correlation
[row,col] = find(C==minC,1,'first')

 %  -0.9441 row = 301 col = 20 :
 
gene_k = gene_names(301,1)%'B0303.2'
gene_l = gene_names(20,1) %'sqt-1'
 
%show appropriate plots of the gene expression profiles (across all time points) for these two pairs of genes. 

x= B(428,:);
y= B(406,:);
figure; scatter(x,y)
title('Positively Correlated Pair of Dynamic Genes'); %title plot
xlabel('csn-2')
ylabel('sra-35')


x= B(301,:);
y= B(20,:);
figure; scatter(x,y)
title('Negatively Correlated Pair of Dynamic Genes'); %title plot
xlabel('B0303.2')
ylabel('sqt-1')



%% Question B) Now, calculate pairwise correlations of all the samples (time points) in the C.elegans dataset

A= C_elegans_RNA_Seq_Embryo;
% Calculate pairwise correlations of these dynamic genes
B= log10(1+A);
% gene_i = B(500,CE_genes);
x = B();
y = B();
C = corr(x',y');
C(find(isnan(C))) = 0;
[~, xi] = sort(C,'descend');
D = xi(1:3394)
% C(xi(1:10)) % find correlation      

C(C==1)=0;

% imagesc(C);

%identify the most positively correlated pair and most negatively correlated pair in this gene set. 
maxC = max(max(triu(C,1))) % highest pos correlation
[row,col] = find(C==maxC,1,'first')

%     0.9989 row = 2238 , col = 6992 :

gene_i = CE_genes(2238,1) %'C18H7.13'
gene_j = CE_genes(6992,1) %'srj-21'



minC = min(min(triu(C,1))) % highest neg correlation
[row,col] = find(C==minC,1,'first')

 %  -0.9656 row = 461 col = 19 :
 
gene_k = CE_genes(461,1)%'baf-1'
gene_l = CE_genes(19,1) %'soc-2'
 
%show appropriate plots of the gene expression profiles (across all time points) for these two pairs of genes. 

x= B(2238,:);
y= B(6992,:);
figure; scatter(x,y)
title('Positively Correlated Pair of All Genes'); %title plot
xlabel('C18H7.13')
ylabel('srj-21')


x= B(461,:);
y= B(19,:);
figure; scatter(x,y)
title('Negatively Correlated Pair of All Genes'); %title plot
xlabel('baf-1')
ylabel('soc-2')

%%

lili_vector = [zeros(1,24), ones(1,10), zeros(1,16)];
    faye_vector = [zeros(1,10), ones(1,10),zeros(1,10), ones(1,10), zeros(1,10)];
        figure; imagesc(faye_vector)
    y = B;
    C = corr(faye_vector',y');
    C(find(isnan(C))) = 0;
    [~, xi] = sort(C,'descend');
    figure; plot(B(xi(2:21),:)')
    figure; imagesc(B(xi(2:200),:))

  y =B;
C = corr(faye_vector',y');
C(find(isnan(C))) = 0;


[~, xi] =sort(C, 'descend'); % find indices
D = xi(1:10)

figure; plot(B(xi(1:10),:))

figure; imagesc(B(xi(1:10),:))  

%% 
%a) Identify 500 dynamic genes in the pancreas single cell dataset.
load('pancreas_data')
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


condition_3= find(sum(B')>28);

dyn_genes= intersect(condition_3,intersect(condition_1,condition_2));length(dyn_genes)

figure;imagesc(A_log(dyn_genes,:))


C=1-squareform(pdist(A_log(dyn_genes,:)', 'correlation')); % added correlation and 1- to find correaltion not euclidian

figure;imagesc(C);

%% Perform a tSNE analysis on the dynamic gene set and present a scatterplot of tSNE1 versus tSNE2.
A_for_pca = log10(1+A(dyn_genes,:));
[COEFF, A_SCORE, LATENT, TSQUARED, EXPLAINED] = pca(A_for_pca');

    rng(1);
    X = A_for_pca';
    labels = [];
    no_dims = 3;
    initial_dims = 4;
    perplexity = 10;
    mappedX = tsne(X, labels, no_dims, initial_dims, perplexity);
    
%     figure;
%     subplot(1,2,1); scatter(mappedX(:,1),mappedX(:,2),10);
%     subplot(1,2,2); scatter(A_SCORE(:,1),A_SCORE(:,2),10);

[~, xi ] = sort(dyn_genes');
%display the gene names of the genes with the highest coefficients
gene_names=gene_names(xi(1:500))


figure;scatter(mappedX(:,1),mappedX(:,2),'filled');
xlabel('tSNE1');
ylabel('tSNE2');
title('tSNE1 vs. tSNE2 by dyn genes');
colorbar;

%% %b) Cluster the pancreas single cell data using k?means clustering,
% change order according to k
km_param=50;
km=kmeans(A_log(dyn_genes,:)',km_param);

% change order according to k

[i,xi]=sort(km);
% figure;imagesc(C)
i = color

C= corr(A_log(dyn_genes,:));
figure;imagesc(C(xi,xi)); %figure;imagesc(C(xi,xi), [0.2 0.5]);


figure;scatter(mappedX(:,1),mappedX(:,2),km,'filled');
xlabel('tSNE1');
ylabel('tSNE2');
title('tSNE1 vs. tSNE2 by kmeans');
colorbar;



