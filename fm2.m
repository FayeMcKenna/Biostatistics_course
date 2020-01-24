% % Homework 3 Solutions
% % Question 1. load('C_elegans_dataset.mat'); load CE_genes.mat
% %part a.
% %find 500 dynamic genes
% C= C_elegans_RNA_Seq_Embryo;
% std_thresh= find(std(C')> 2); %want standard deviation larger than 2 sum_thresh= find(sum(C')>5000); %sum of expression of genes >5000 combo= intersect(sum_thresh, std_thresh);
% D= log10(1+C); %log scaled vector of original A values
% dyn_genes= D(combo,:);
% [s, i]=sort(sum(dyn_genes'), 'descend');
% dyn_genes_final= dyn_genes(i(1:500),:);
% % correlations
% corrs= corr(dyn_genes_final'); corrs= tril(corrs, -1);
% % find most positively and negatively correlated [max_corr, max_ind]= max(corrs(:));
% [min_corr, min_ind]= min(corrs(:));
% [max_row, max_col] = ind2sub(size(corrs),max_ind); [min_row, min_col] = ind2sub(size(corrs),min_ind);
% names_all= CE_genes;
% dynamic_names= names_all(i(1:500));
% max_gene1= dynamic_names(max_row);
% max_gene2= dynamic_names(max_col);
% min_gene1= dynamic_names(min_row);
% min_gene2= dynamic_names(min_col); names=[max_gene1, max_gene2, min_gene1, min_gene2];
% figure;
% subplot(1,2,1);
% scatter(dyn_genes_final(max_row,:), dyn_genes_final(max_col,:), 10, 'k', 'filled'); xlabel('B0379.7 Expression');
% ylabel('B0041.11 Expression');
% title('Most Positive Correlation Between 2 Dynamic Genes');
% subplot(1,2,2);
% scatter(dyn_genes_final(min_row,:), dyn_genes_final(min_col,:), 10, 'k', 'filled'); xlabel('B03261.5 Expression');
% ylabel('C07G3 Expression');
% title('Most Negative Correlation Between 2 Dynamic Genes');

 Y = xlsread('Y.xlsx')
outpath = '/Users/fayemckenna/Desktop/Sackler/semester2/biostats/EDA';
outname = 'fmfig';

heatscatter(X, Y, outpath, outname)
%heatscatter(X, Y, outpath, outname, numbins, markersize, marker, plot_colorbar, plot_lsf, xlab, ylab, title)
