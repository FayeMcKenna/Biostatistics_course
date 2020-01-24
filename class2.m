% PCA- reduces dimensionality 

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
    load pancreas_data.mat
end    

% figure; imagesc(log10(1+A)) % put to scale in matrix
% 
% %pca can use explained to see how much of the variation was explained
% 
% A_for_pca=log10(1+A);
% [COEFF, SCORE]=pca(A_for_pca);
% 
% height =[72 75 60 80 55 48];
% weight = [78 70 60 85 60 40];
% %figure; plot(height, weight, 'ko');
% 
% % make one dimenstion of the two variables
% B= [height; weight];
% [B_COEFF, B_SCORE]= pca(B');

% B score will be ranked by more effect 
% what you give it you get back but then you decide what you care about 

% second fun through

% how much of variance is explained? may need to play around with transpose

A_for_pca=log10(1+A);
[A_COEFF, A_SCORE,latent,tsquared,explained,mu]=pca(A_for_pca');

% now the output is like a new gene represented of pca
figure; plot(A_SCORE(:,1),A_SCORE(:,2),'.')

% x is pca1, y is pca2: look at coefficients which genes are contributing
% to the axis


% how much of variance is explained? numbers go down in a row
explained(1:2)

% can color dots by influence of insulin

insulin_gene=strmatch('INS',gene_names,'exact');

%figure; plot(A_SCORE(:,1), A_SCORE(:,2), '.')
% scatter you can control size and color of dots

figure; scatter(A_SCORE(:,1), A_SCORE(:,2),5,A_for_pca(insulin_gene,:), 'filled')
% 5 is size color is expression level

%subplot

%othergenes
gcg_gene=strmatch('GCG',gene_names,'exact');
ppy_gene=strmatch('PPY',gene_names,'exact');
sst_gene=strmatch('SST',gene_names,'exact');

figure
scatter(A_SCORE(:,1), A_SCORE(:,2), 5, A_for_pca(gcg_gene,:), 'filled')

figure;
subplot(2,2,1)
scatter(A_SCORE(:,1), A_SCORE(:,2),5,A_for_pca(insulin_gene,:), 'filled')
% 5 is size color is expression level


subplot(2,2,2)
scatter(A_SCORE(:,1), A_SCORE(:,2),5,A_for_pca(ppy_gene,:), 'filled')
% 5 is size color is expression level


subplot(2,2,3)
scatter(A_SCORE(:,1), A_SCORE(:,2),5,A_for_pca(sst_gene,:), 'filled')
% 5 is size color is expression level


subplot(2,2,4)
scatter(A_SCORE(:,1), A_SCORE(:,2),5,A_for_pca(gcg_gene,:), 'filled')
% 5 is size color is expression level


% first 10 can explain how much?
sum(explained(1:10))

% size coefficients

size(A_COEFF)
figure;hist(A_COEFF(:,1),1000) % 1000 bins
% coeffecicients in pca: what are they? if you are a cell you have 1257
% numbers, where is a cell on pC one? look at expression and coefficients
% and sum what you are for pc1
% the only genes that make a difference of genes that are either high or
% low not zero, you can see in histogram that many are multiplied by zero


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
%now we can see that pc1 has to do with beta cells are you or are you not
%beta cells?
% dont care about direction (absolute value) but just that there is a
% relationship, can narrow down on most expressed ex then your pca willbe
% very explained
thresh = 500; std_thresh= 2;
expressed=find(sum(A')>thresh); %works on column so have to transpose
some_variation=find(std(A')>std_thresh);
dyn_genes= intersect(expressed, some_variation);

A_for_pca=log10(1+A(dyn_genes,:));
[A_COEFF, A_SCORE, latent, tsquared, explained, mu]=pca(A_for_pca');

%othergenes
gcg_gene=strmatch('GCG',gene_names,'exact');
ppy_gene=strmatch('PPY',gene_names,'exact');
sst_gene=strmatch('SST',gene_names,'exact');
insulin_gene=strmatch('INS',gene_names,'exact');


figure
scatter(A_SCORE(:,1), A_SCORE(:,2), 5, A_for_pca(gcg_gene,:), 'filled')

% figure;
subplot(2,2,1)
scatter(A_SCORE(:,1), A_SCORE(:,2),5,A_for_pca(insulin_gene,:), 'filled')
% 5 is size color is expression level


subplot(2,2,2)
scatter(A_SCORE(:,1), A_SCORE(:,2),5,A_for_pca(ppy_gene,:), 'filled')
% 5 is size color is expression level


subplot(2,2,3)
scatter(A_SCORE(:,1), A_SCORE(:,2),5,A_for_pca(sst_gene,:), 'filled')
% 5 is size color is expression level


subplot(2,2,4)
scatter(A_SCORE(:,1), A_SCORE(:,2),5,A_for_pca(gcg_gene,:), 'filled')
% 5 is size color is expression level
% 
% 
%% class 3 continued

% how do you know transpose or not with pca?

% with pca, have variables and cell ; for example genes (rows) are
% variables but they want the rows to be observations and columns cell are
% observations

% review what does explain do?

explained(1:10) % 14% of variation is explained by first pca 
% if you want to look at more dimensions use tesne, it is possible to do 3d

%% class 3

% part_1_tsne
% A_for_pca=log10(1+A(dyn_genes,:));
% [A_COEFF, A_SCORE, LATENT, TSQUARED, EXPLAINED]=pca(A_for_pca');
% 
% X= A_for_pca';
% labels=[];
% no_dims=3;
% initial_dims = 10;
% perplexity= 30;
% mappedX= tsne(X, labels, no_dims, initial_dims, perplexity);
% 
% figure; 
% subplot(1,2,1);scatter(mappedX(:,1),mappedX(:,2),10);
% subplot(1,2,2);scatter(A_SCORE(:,1),A_SCORE(:,2),10); 

% wont be the same everytime
% perplexity is the one tuning variable

%save('mappedx.mat') ; % save because everytime something new 
% could also use rng(1);above the command and will run same thing everytime

%part2_
load('C_elegans_dataset.mat')
load('CE_genes.mat')

B= log10(1+C_elegans_RNA_Seq_Embryo); 
gene_i = strmatch('dsl-1', CE_genes);
figure; plot(B(gene_i,:));

% how do you find the genes that have the same expression pattern as the
% plot? Euclidian distance and correlation
%correlation
x=[1 2 2 3 3 4 5 5.5]; y = [1 .8 2 2.3 4 3.7 5.5 4.2;];

%find correlation by % really just standardize x and standardize y mult
%together and div length 
corr(x',y')

% now find other profile that is same by correlating it with all other
% profiles of genes

gene_j=strmatch('mex-3', CE_genes); % expression profile of these two genes

x= B(gene_i,:);
y= B(gene_j,:);
figure; scatter(x,y)

% do them all 
x =B(gene_i,:);
y = B;
C=corr(x',y')

C(find(isnan(C))) = 0;


[~, xi] =sort(C, 'descend'); % find indices
D = xi(1:10)

C(xi(1:10)) % find correlation

gene_j= 9857; % expression profile of this high corr gene

x= B(gene_i,:);
y= B(gene_j,:);
figure; scatter(x,y)

% you lose the direction when you sort which you may want


% now find opposite
xi(end-10:end)

% look at ones all most correlated
gene_j= 9857; % expression profile of this high corr gene

x= B(gene_i,:);
y= B(gene_j,:);
figure; scatter(x,y)

% [~, xi] = sort (C, 'descend'];
figure; plot(B(xi(2:21),:)')

% how do you find the most correlations , makeyour own vector

lili=[zeros(1,24), ones(1,10), zeros(1,16)];
figure; imagesc(lili)

% our own vector
y =B;
C = corr(lili',y');
C(find(isnan(C))) = 0;


[~, xi] =sort(C, 'descend'); % find indices
D = xi(1:10)

figure; plot(B(xi(2:21),:))

figure; imagesc(B(xi(2:21),:))

% compare both created vectors

%% euclidian distance

gene_j=strmatch('mex-3', CE_genes); % expression profile of these two genes

x= B(gene_i,:);
y= B(gene_j,:);
figure; scatter(x,y)

%pdist

for i=1:20517
    y = B(i,:);
    D(i)=pdist([x;y]); % save distance
end

[~,xI] = sort(D);
figure;plot(B(xi(2:21),:)')
figure; imagesc(B(xi(2:200)))


% no nans because you are never dividing 
% euclidian is you want low, you want closest 
% sort by accend 
% pearsons is standardize so it looks less tight, and pearsons has a range
% which gives you more information 




