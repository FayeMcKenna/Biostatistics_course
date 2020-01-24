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
    load cortex_data.mat
end    

% Question 1: A gene is detected in a cell if its expression is >0.

%1A) Calculate the number of genes detected for each cell

A_log=log10(1+A); % log of data 
% figure; imagesc(A_log);

B= A_log; B(find(B>0))=1; % find where gene expression in data set is greater than zero, set that to 1

for j =1:1000 % for each column or 'cell'
    a=B(:,j);
    sumGene(j) = sum(a);  % count each '1' of gene expressionin column or 'cell'
end

% 1B) Calculate the mean, median and standard deviation for the sum of
% all genes across all cells 

% mean, med, stdv of sum of genes
 meanGene = mean(sumGene); % mean of of sum of all genes
 medGene = median(sumGene); % median of of sum of all genes
 stdGene = std(sumGene); % standard deviation of sum of all genes
   
%1C) Create a histogram of the number of genes detected for all cells with
%appropriate num of bins

% find number of genes detected in all cells
for k =1:19972 % for each gene
    a=B(k,:); % in row of B matrix
    E(k) = sum(a); % sum expression of genes in each cell
end

% filter to take out zeros from E
gEnes = E;
b=size(gEnes);
gEnes(gEnes==0) = [];

% plot histogram of gene expression 
figure; hist(log10(1+gEnes(1,:)),40);%scale 50 bins
xlabel('Log10 Gene Expression Levels for all Cells');
ylabel('Freq of Genes in cells');
title('Number of Genes Detected for all Cells');

% Question 2: Choose the top 500 dynamic genes and explain how you picked the
% dynamic genes and how you set the thresholds

% I selected the dynamic genes that were expressed, highly expressed, and had some
% variation (std threshold).

% condition_1=expressed (give top half)
condition_1=find(sum(A')>median(sum(A'))); length(condition_1)

% condition_2=standard deviation greater than 3
condition_2=find(std(A')>3);

% condition_3= number of cells expressed in , take the matrix whatever is
% higher than one sum and greater than threshold
B= A_log; B(find(B>0))=1;

condition_3= find(sum(B')>701); 

% dynamic genes have to meet all three conditions
dyn_genes= intersect(condition_3,intersect(condition_1,condition_2));length(dyn_genes)

%plot in heatmap the 500 dynamic genes unsorted
figure;imagesc(A_log(dyn_genes,:))
title('Heatmap of 500 most dynamic genes'); %title plot
xlabel('cells')
ylabel('Dyn Genes')
colorbar; %include colorbar

y =(A_log(dyn_genes,:));
%plot in heatmap the 500 dynamic genes sorted on variable cluster ID

% [~, xi ]=sort(A', 'clust_ID');
C= [clust_ID,(A)']; % add two matrixes together in order to sort
[~,idx] = sort(C(:,1)); % sort just the first column of combined matrix 
sortedmat = C(idx,:);   % sort the whole matrix using the first cluster_ID sorted indices

D = transpose(sortedmat(:,2:end)); % recreate matrix in order of cluster_ID sorting

D_log=log10(1+D); % log of new matrix data 

figure;imagesc(D_log(dyn_genes,:))
title('Heatmap of 500 most dynamic genes in cluster_ID order'); %title plot
xlabel('cells')
ylabel('Dyn Genes')
colorbar; %include colorbar



% the difference between the two is that the first is ordered based on the
% cells in a understandibly random order, the second is ordered based on
% low to high clustering ID

%Question 3: Perfor `m PCA. Remember to scale the data and to lookonly at the dynamic genes

% % create new matrix of 
% G= transpose(1:19972);
% G = [G,A]; % new matrix with order of rows
% dyn_genes=dyn_genes'; % transpose gyn_genes rows
% 
% [a_sorted, a_order] = sort(dyn_genes); % find rows that are that of dyn_genes
% newG = G(dyn_genes,:);
% 
% G=newG(:,2:end); % remove row label from dyn_genes

% create new matrix of 
G= transpose(1:19972);
G = [G,A]; % new matrix with order of rows
dyn_genes=dyn_genes'; % transpose gyn_genes rows % transpose gyn_genes rows

[a_sorted, a_order] = sort(dyn_genes); % find rows that are that of dyn_genes
newG = G(dyn_genes,:);

y=newG(:,2:end); % remove row label from dyn_genes

%new list of gene names

% m=dyn_genes;
% 
% [I]= cell(gene_names);
% 
% for j=m(1,:)
%     [new_gn]=cell(I(j)) 
% end

%compute PCA

 A_for_pca=log10(1+y);
[A_COEFF, A_SCORE,latent,tsquared,explained,mu]=pca(A_for_pca);

%3A)Figure of PCA 1 vs PCA 2
figure; plot(A_SCORE(:,1),A_SCORE(:,2),'.')
title('PCA1 vs PCA2 of dynamic genes'); %title plot
xlabel('PCA1')
ylabel('PCA2')

explained(1:2)

%3B) Color the scatter plot from part (a) based on the following 4 marker genes: ?Gad1?,?Tbr1?,?Spink8? and ?Mbp?
%display a plot for each gene (use the subplot function).


A_for_pca=log10(1+A);
[A_COEFF, A_SCORE,latent,tsquared,explained,mu]=pca(A_for_pca');
% A_for_pca=transpose(A_for_pca);
figure; plot(A_SCORE(:,1),A_SCORE(:,2),'.')
title('PCA1 vs PCA2 of dynamic genes'); %title plot
xlabel('PCA1')
ylabel('PCA2')



Gad_gene=strmatch('Gad1',gene_names,'exact');
Tbr_gene=strmatch('Tbr1',gene_names,'exact');
Spink_gene=strmatch('Spink8',gene_names,'exact');
Mbp_gene=strmatch('Mbp',gene_names,'exact');

figure; plot(A_SCORE(:,1),A_SCORE(:,2),'.')
title('PCA1 vs PCA2 of dynamic genes'); %title plot
xlabel('PCA1')
ylabel('PCA2')

tApca=A_for_pca;
%tApca=transpose(A_for_pca);

%plot 4 genes with subplots
figure;
subplot(2,2,1)
scatter(A_SCORE(:,1), A_SCORE(:,2),5,tApca(Mbp_gene,:), 'filled')
% 5 is size color is expression level


subplot(2,2,2)
scatter(A_SCORE(:,1), A_SCORE(:,2),5,tApca(Spink_gene,:), 'filled')
% 5 is size color is expression level


subplot(2,2,3)
scatter(A_SCORE(:,1), A_SCORE(:,2),5,tApca(Tbr_gene,:), 'filled')
% 5 is size color is expression level

subplot(2,2,4)
scatter(A_SCORE(:,1), A_SCORE(:,2),5,tApca(Mbp_gene,:), 'filled')

% use Clus ID to code map
figure
subplot(1,2,1);scatter(A_SCORE(:,1), A_SCORE(:,2),10,sortedmat(:,1), 'filled')


%3D)

%PC1= ??
PC1=A_COEFF(:,1);

%absolute value
aPC1=abs(PC1)
%sort them and remember indexes
[~, xi ]=sort(aPC1, 'descend'); % sorting absolute with highest first and getting the index
%display the gene names with the highest coefficients
highest_gene=gene_names(xi(1:10)) 
%now we can see that pc1 has to do with beta cells are you or are you not
%beta cells?
%look at PC2 psuedo code
PC1=A_COEFF(:,2);

%absolute value
aPC1=abs(PC1)
%sort them and remember indexes
[~, xi ]=sort(aPC1, 'descend'); % sorting absolute with highest first and getting the index
%display the gene names with the highest coefficients
highest_gene=gene_names(xi(1:10)) 

%PLp1

%absolute value
aPC1=abs(PC1)
%sort them and remember indexes
[~, xi ]=sort(aPC1, 'ascend'); % sorting absolute with highest first and getting the index
%display the gene names with the highest coefficients
highest_gene=gene_names(xi(1:10)) 
%now we can see that pc1 has to do with beta cells are you or are you not
%beta cells?
%look at PC2 psuedo code
PC1=A_COEFF(:,2);

%absolute value
aPC1=abs(PC1)
%sort them and remember indexes
[~, xi ]=sort(aPC1, 'ascend'); % sorting absolute with highest first and getting the index
%display the gene names with the highest coefficients
highest_gene=gene_names(xi(1:10)) 

%'Fabp4'


plp_gene=strmatch('Plp1',gene_names,'exact');
fab_gene=strmatch('Fabp4',gene_names,'exact');


subplot(2,2,1)
scatter(A_SCORE(:,1), A_SCORE(:,2),5,tApca(plp_gene,:), 'filled')
% 5 is size color is expression level

subplot(2,2,2)
scatter(A_SCORE(:,1), A_SCORE(:,2),5,tApca(fab_gene,:), 'filled')


%Question 4 % do TSNE

    rng(1);
    X = A_for_pca';
    labels = [];
    no_dims = 3;
    initial_dims = 4;
    perplexity = 10;
    mappedX = tsne(X, labels, no_dims, initial_dims, perplexity);
    
    figure;
    plot(( scatter(mappedX(:,1),mappedX(:,2),10,clust_ID, 'filled'))
    
     subplot(1,2,1); scatter(mappedX(:,1),mappedX(:,2),10,sortedmat(:,1));
%     subplot(1,2,2); scatter(A_SCORE(:,1),A_SCORE(:,2),10);


