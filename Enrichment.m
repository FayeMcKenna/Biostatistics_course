function [pVals, H,ints] = Enrichment(gene_list, p_thresh)
import bioma.*;    import bioma.data.*;    import statistics.*

load Gene_ontology_matrix.mat;    A = Gene_ontology; names = Gene_ontology_names;
G = zeros( length(A) , 1);
gene_list
G(gene_list,1) = 1;

all_genes = length(A);          
for j = 1:size(G,2)    
    for k = 1:size(A,2)
        a = length(intersect(find(G(:,j)),find(A(:,k)))); b = all_genes; c = sum(G(:,j)); d = sum(A(:,k));
        ints(k) = a;
        pVals(j,k) = hygecdf(c-a-1,b,c,b-d);
    end
end

i=find(pVals < p_thresh);
H = DataMatrix(pVals(:,i),{''},names(i));