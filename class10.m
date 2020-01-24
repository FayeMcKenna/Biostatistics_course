threshold = 2;

try
    CE_genes;
catch
    load CE_genes.mat
    load C_elegans_dataset.mat
end

%% Regression
rng(4);
samples_to_exclude = randi(50,1,2);

A = log10(1+C_elegans_RNA_Seq_Embryo(:, setdiff(1:50, samples_to_exclude)));

highly_expressed = find(max(A')>threshold);
B = A(highly_expressed,:);

%make the profile
profile = 1:48;

%correlate said profile with each gene
R = corr(profile',B');

%sort so as to identify the gene that best correlates with the profile

[i,xi] = max(R)
% figure; plot(B(xi,:),'ro');

X = setdiff(1:50, samples_to_exclude);
Y = B(xi,:);
model = fitlm(Y, X);
% figure; plot(model)

values = log10(1+C_elegans_RNA_Seq_Embryo(highly_expressed(xi),samples_to_exclude));
intercept = -7.9665; x1 =  13.14 ;
intercept + values*x1
samples_to_exclude



