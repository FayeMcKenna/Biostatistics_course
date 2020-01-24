normalization_part = 0;
ductal_cell_part   = 1;

if ductal_cell_part
   try
       A_ductal;
   catch
       load working_data_ductal.mat;
   end

   high_mean = find(mean(A_ductal')>50);
   high_var = find(std(A_ductal')>1000);
   dyn_exp = intersect(high_mean, high_var);
   length(dyn_exp)
   B = log10(1+A_ductal);
%     figure; imagesc(B);    

   [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(B(dyn_exp,:)');
   figure; scatter(SCORE(:,1),SCORE(:,2),'o');
   
   group1 = find(SCORE(:,1)>-1)
   group2 = find(SCORE(:,1)<=-1);
   
   figure;
   imagesc(B(dyn_exp,[group1' group2']))
   
   g1 = dyn_exp(1);
   [~, p] = ttest2(B(g1,group1),B(g1,group2));
   
   [i,xi] = sort(abs(COEFF(:,1)));
   gene_names(dyn_exp(xi))
end

if (normalization_part)
   try
       CE_genes;
   catch
       load C_elegans_dataset.mat
       load C_elegans_dataset_raw.mat
       load CE_genes.mat
   end
   
   %  figure; hist(sum(C_elegans_RNA_Seq_Embryo_Raw),100)
   rng(1); M = rand(10,5);
   M_tpm = [M./mean(M)].*1000000;
   
   %if there is time discuss quantile normalization
   
   A_tpm = [C_elegans_RNA_Seq_Embryo_Raw./sum(C_elegans_RNA_Seq_Embryo_Raw)].*1000000;
   isequal(C_elegans_RNA_Seq_Embryo, A_tpm);
   
   
   %identify expressed genes
   highly_expressed = find(max(A_tpm')>100);
   have_some_lows  = find(min(A_tpm')==0);
   dyn_exp = intersect(highly_expressed,have_some_lows );
   
   %sort them
   [~, T] = max(A_tpm(dyn_exp,:)');
   [i,xi] = sort(T);
   
   %view it
   figure;
   subplot(1,2,1);
   imagesc(log10(1+A_tpm(dyn_exp(xi),:)));
   
   %z-score it and view alongside it
   
   B = log10(1+A_tpm(dyn_exp(xi),:));
   B_z = zscore(B')';
   subplot(1,2,2);
   imagesc(B_z,[0 2]);
   
   
   % v = rand(1,10);
   % v_z = [v-mean(v)]./std(v)
   
   
   %talk in a second about why it's called a Z-score
end


