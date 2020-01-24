part_1 = 0;
part_2 = 0;
part_3 = 1;


if (part_3)
    
    P_val_ref = [0.05*[1:532]]/532;
    [i,xi] = sort(P_vals);
    P_vals_sorted = P_vals(xi);

    figure; plot(P_val_ref); hold on;
    plot(P_vals_sorted); hold on;
    
    FDR(P_vals, 0.05)
    length(find(P_vals<FDR(P_vals, 0.05)))
    length(find(P_vals<0.05))

end

if (part_2)
    %define a profile
    
    beck_prof = [zeros(1,10), ones(1,10), zeros(1,10), ones(1,10), zeros(1,10)];
    itai_prof = [zeros(1,5), ones(1,5), zeros(1,40)];
    
    %decide on the highly expressed
    load C_elegans_dataset.mat
    load CE_genes.mat
    
    A = log10(1+C_elegans_RNA_Seq_Embryo);
    B = zscore(A')';
    
    highly_expressed = find(mean(C_elegans_RNA_Seq_Embryo')>500);
    
    %correlate
    C = corr(itai_prof',C_elegans_RNA_Seq_Embryo(highly_expressed,:)');
    
    %sort
    [i,xi] = sort(C','descend');
    
    %pick the most correlated
    gene_list = highly_expressed(xi(1:30));
    
    %compute the enrichment
    [P_vals,a] = Enrichment(gene_list,0.0000001)
    
    
end


if (part_1)
    
    %simulating data
    
    groupings = [ones(1,50),zeros(1,50)];
    clear P;
    for i = 1:1000
        measurements = randn(1,100);
        [h,P(i)] = ttest2(measurements(find(groupings==1)),measurements(find(groupings==0)));
    end
    length(find(P<0.05))
end

