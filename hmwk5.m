%% clean up %Hmk3
clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same

load C_elegans_dataset.mat
load('CE_genes.mat')
load('Gene_ontology_matrix.mat')


A= C_elegans_RNA_Seq_Embryo;
%Calculate number of genes for each cell
A_log= log10(1+A);
index=A_log~=0;
Num_of_genes=sum(index,1);

A = log10(1+C_elegans_RNA_Seq_Embryo);

%% part 1 and 2
%vector pattern A
a    = zeros(1,50); %patterna
%beggining low
for i=1:16 
    a(i)  = 10;                     
end
%middle high
for i=17:35 
    a(i)    = 45;                       
end
%end low
for i = 36:50 
    a(i)=10;
end

%pattern A : low beginning, high middle, low end
patternA = a';
corr_patternA= corr(A_log',patternA);
% find 30 genes most correlated to this pattern
corr_patternA(find(isnan(corr_patternA))) = 0;
[corr_pattern_sortedA,i] = sort(corr_patternA,'descend');
max_patternA_idx = i(1:30);
max_patternA_names = CE_genes(max_patternA_idx);
figure; subplot(1,2,1); plot(patternA,'Color','k','Linewidth',3); title('Pattern')
subplot(1,2,2); plot(A_log(max_patternA_idx,:)','Linewidth',2); legend(max_patternA_names);
title('Top 30 genes correlated to the pattern A');
xlabel('Time Point');
ylabel('Gene Expression');

%pattern B : low beginning high end
patternB = [1:50]';
corr_patternB= corr(A_log',patternB);
% find 30 genes most correlated to this pattern
corr_patternB(find(isnan(corr_patternB))) = 0;
[corr_pattern_sortedB,j] = sort(corr_patternB,'descend');
max_patternB_idx = j(1:30);
max_patternB_names = CE_genes(max_patternB_idx);
figure; subplot(1,2,1); plot(patternB,'Color','k','Linewidth',3); title('Pattern')
subplot(1,2,2); plot(A_log(max_patternB_idx,:)','Linewidth',2); legend(max_patternB_names);
title('Top 30 genes correlated to the pattern B');
xlabel('Time Point');
ylabel('Gene Expression');

%pattern C : high beginning low end
patternC = flipud(patternB);
corr_patternC= corr(A_log',patternC);
% find 30 genes most correlated to this pattern
corr_patternC(find(isnan(corr_patternC))) = 0;
[corr_pattern_sortedC,k] = sort(corr_patternC,'descend');
max_patternC_idx = k(1:30);
max_patternC_names = CE_genes(max_patternC_idx);
figure; subplot(1,2,1); plot(patternC,'Color','k','Linewidth',3); title('Pattern')
subplot(1,2,2); plot(A_log(max_patternC_idx,:)','Linewidth',2); legend(max_patternC_names);
title('Top 30 genes correlated to the pattern C');
xlabel('Time Point');
ylabel('Gene Expression');

%% part 3 Gene enrichment analysis
%simulating data, calculate using alpha bonferroni

% bonferroni P calculation
Bo =.01/532; 
   
% Define genes of pattern A
gene_listA = max_patternA_idx';  % indices of gene list
    
%compute the enrichment
[P_valsA,a] = Enrichment(gene_listA,Bo) % P calculated
   
% FDR correction
P_val_refA = [.01 *[1:532]]/532;
[i,xi] = sort(P_valsA);
P_vals_sortedA = P_valsA(xi);

figure; plot(P_val_refA); hold on;
plot(P_vals_sortedA); hold on;
    
FDR(P_valsA, .01)
length(find(P_valsA<FDR(P_valsA, .01)))
length(find(P_valsA<.01))

%compute the enrichment FDR
[P_valsA2,a2] = Enrichment(gene_listA,4.6884e-04) % P calculated

%% B
%Define genes of pattern B
gene_listB = max_patternB_idx';  % indices of gene list
    
%compute the enrichment
[P_valsB,b] = Enrichment(gene_listB,Bo) % P caluculated
   
% FDR correction
P_val_refB = [0.01*[1:532]]/532;
[i,xi] = sort(P_valsB);
P_vals_sortedB = P_valsB(xi);

figure; plot(P_val_refB); hold on;
plot(P_vals_sortedB); hold on;
    
FDR(P_valsB, 0.01)
length(find(P_valsB<FDR(P_valsB, 0.01)))
length(find(P_valsB<0.01))

%compute the enrichment FDR
[P_valsB2,b2] = Enrichment(gene_listB,3.1904e-04) % P caluculated

%% C
%Define genes of pattern C
gene_listC = max_patternC_idx';  % indices of gene list
    
%compute the enrichment
[P_valsC,c] = Enrichment(gene_listC,Bo) % P calculated
   
% FDR correction
P_val_refC = [0.01*[1:532]]/532;
[i,xi] = sort(P_valsC);
P_vals_sortedC = P_valsC(xi);

figure; plot(P_val_refC); hold on;
plot(P_vals_sortedC); hold on;
    
FDR(P_valsC, 0.01)
length(find(P_valsC<FDR(P_valsC, 0.01)))
length(find(P_valsC<0.01))

%compute the enrichment FDR
[P_valsC2,c2] = Enrichment(gene_listC,8.1879e-04) % P calculated

%% part 3 Visualize data 

% bonferonni organize data for plotting
% pattern A genes and Go categories
Ax = double(a)'
Ay = [zeros(1,23)]';
Ax = [Ax;Ay]

A_names = colnames(a)'

% pattern B genes and Go categories
Bx = double(b)'
Bz = [zeros(1,2)]';
Bk = [zeros(1,21)]';
Bx = [Bz;Bx;Bk]

B_names = colnames(b)'

% pattern C genes and Go categories
Cx = double(c)'
Cy = [zeros(1,4)]'
Cx = [Cy;Cx]

C_names = colnames(c)'

Bon_name = [ A_names; B_names; C_names]

save('Bon_name')

bar(Ax) 

x = categorical({names1}); y = [alpha1]; b = bar(x,y); 
title('Breakdown of the average pump power consumption') 
ylabel('Average Power (kW)') 
xlabel('Month') 
labels = {'Overall','Fey Pivot','Pods'}; 
legend(labels,'Location','southoutside','Orientation','horizontal')



%% FDR organize data for plotting
% pattern A genes and Go categories
Ax = double(a2)'
Ay = [zeros(1,61)]';
Ax = [Ax;Ay]

A_names = colnames(a2)'

% pattern B genes and Go categories
Bx = double(b2)'
Bz = [zeros(1,25)]';
Bk = [zeros(1,41)]';
Bx = [Bz;Bx;Bk]

B_names = colnames(b2)'

% pattern C genes and Go categories
Cx = double(c2)'
Cy = [zeros(1,45)]'
Cx = [Cy;Cx]

C_names = colnames(c2)'

Bon_name = [ A_names; B_names; C_names]

bar(Ax) 

x = categorical({names1}); y = [alpha1]; b = bar(x,y); 
title('Breakdown of the average pump power consumption') 
ylabel('Average Power (kW)') 
xlabel('Month') 
labels = {'Overall','Fey Pivot','Pods'}; 
legend(labels,'Location','southoutside','Orientation','horizontal')

