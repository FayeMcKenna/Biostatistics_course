%% clean up
clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same

%Question1
% Create magic matrix 50x50
 magMat=magic(50); 

 % Create heatmap of magic matrix
figure; imagesc(magMat);

%Question2
% find mean, max, min, median, and standard deviation each row in magic
% matrix

for j =1:50
    a=magMat(j,:)
    meanRow(j) = mean(a)
    maxRow(j) = max(a);
    minRow(j) = min(a);
    medRow(j) = median(a);
    stdRow(j) = std(a);
    outputrow = [meanRow; maxRow; minRow; medRow; stdRow]      
end
  

% find mean, max, min, median, and standard deviation of each column in magic
% matrix

for j =1:50
    a=magMat(j,:)
    meanRow2(j) = mean(a)
    maxRow2(j) = max(a);
    minRow2(j) = min(a);
    medRow2(j) = median(a);
    stdRow2(j) = std(a);
    outputcol = [meanRow2; maxRow2; minRow2; medRow2; stdRow2]      
end
  
% find mean, max, min, median, and standard deviation of entire magic
% matrix

 meanMat = mean(magMat)
 maxMat = max(magMat);
 minMat = min(magMat);
 medMat = median(magMat);
 stdMat = std(magMat);
 outputMat = [meanMat; maxMat; minMat; medMat; stdMat]    
 
 % present data
%Column Matrix data
 
[~, stage] = max(outputcol'); % find columns index of max value
 
[i,order]= sort(stage); % find indices of sorted 
 
figure; imagesc(outputcol(order,:));


 set(gca, 'ytick', 1:length(outputcol));
 set(gca, 'yticklabel', outputcol(order));
 set(gca, 'yticklabel', {'mean' ,'std' ,'med', 'min' ,'max'});
 
% row matrix data

toutputrow=transpose(outputrow)

[~, stage] = max(toutputrow'); % find columns index of max value
 
[i,order]= sort(stage); % find indices of sorted 
 
figure; imagesc(toutputrow(order,:));
 set(gca, 'xtick', 1:length(outputcol));
 set(gca, 'xticklabel', toutputcol(order));
 set(gca, 'xticklabel', {'mean' ,'std' ,'med', 'min' ,'max'});
 
 % all matrix data
 
 [~, stage] = max(outputcol'); % find columns index of max value
 
[i,order]= sort(stage); % find indices of sorted 
 
figure; imagesc(outputMat(order,:));

 set(gca, 'ytick', 1:length(outputMat));
 set(gca, 'yticklabel', outputMat(order));
 set(gca, 'yticklabel', {'mean' ,'std' ,'med', 'min' ,'max'});
 
% row matrix data

% show histogram of one sample in the C_elegans dataset

% select first col in dataset as sample
data1=C_elegans_RNA_Seq_Embryo(:,1); 

%create list of gene number to plot against
genes=(1:20517)

%plot genes of sample
figure
plot(genes,data1)

% pick 5 genes from CE_genes mat and with high expression over 200 and box

%set array of each gene of interest
gene1=C_elegans_RNA_Seq_Embryo(127,:)
gene2=C_elegans_RNA_Seq_Embryo(280,:)
gene3=C_elegans_RNA_Seq_Embryo(190,:)
gene4=C_elegans_RNA_Seq_Embryo(459,:)
gene5=C_elegans_RNA_Seq_Embryo(18641,:)
genesMat = [gene1; gene2; gene3; gene4; gene5]  

%transpose to plot
tgenesMat= transpose(genesMat)

% plot with boxplot 5 genes
figure
boxplot(tgenesMat)




