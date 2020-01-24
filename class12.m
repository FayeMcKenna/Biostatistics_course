%define a profile

beck_prof = [zeros(1,10), ones(1,10), zeros(1,10), ones(1,10), zeros(1,10)];

%decide on the highly expressed
load C_elegans_dataset.mat
load CE_genes.mat
highly_expressed = find(mean(C_elegans_RNA_Seq_Embryo')>500);

%correlate
C = corr(beck_prof',C_elegans_RNA_Seq_Embryo(highly_expressed,:)');

%sort
[i,xi] = sort(C','descend');

%pick the most correlated
gene_list = highly_expressed(xi(1:30));

%compute the enrichment
Enrichment(gene_list,0.000000001)

if (0)
   (nchoosek(30,6)* nchoosek(70,4)  ) ./  nchoosek(100,10)
   
   
   p6 = (nchoosek(30,6)* nchoosek(70,4)  ) ./  nchoosek(100,10)
   p7 = (nchoosek(30,7)* nchoosek(70,3)  ) ./  nchoosek(100,10)
   p8 = (nchoosek(30,8)* nchoosek(70,2)  ) ./  nchoosek(100,10)
   p9 = (nchoosek(30,9)* nchoosek(70,1)  ) ./  nchoosek(100,10)
   p10 = (nchoosek(30,10)* nchoosek(70,0)  ) ./  nchoosek(100,10)
   
   p6_or_more = p6+p7+p8+p9+p10
   
   
   sum_p = 0;
   for i = 6:10
       p = (nchoosek(30,i)* nchoosek(70,10-i)  ) ./  nchoosek(100,10);
       sum_p = sum_p + p;
   end
   sum_p
   
   
   p = hygepdf(6,100,30,10)
   
   hygepdf(6,100,30,10)
   hygecdf(6,100,30,10)
   hygecdf(0:10,100,30,10)
   figure; plot(hygecdf(0:10,100,30,10))
   
   1 - hygecdf(5,100,30,10)
   
   
   
   
   all_genes = length(A);              % all genes, const
   for j = 1:size(G,2) % do for all custom groups
       for k = 1:size(A,2)
           a = length(intersect(find(G(:,j)),find(A(:,k))));
           b = all_genes;
           c = sum(G(:,j));
           d = sum(A(:,k));
           ints(k) = a;
           pVals(j,k) = hygecdf(c-a-1,b,c,b-d);
       end
   end
end

