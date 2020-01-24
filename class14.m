try
    A;
catch
    load pancreas_data.mat
end

part_1 = 0;
part_2 = 0;
part_3 = 1;

if (part_1)
    marker = strmatch('INS',gene_names,'exact');
    figure; hist(log10(1+A(marker,:)))
    
    x = log10(1+A(marker,:));
    y = zeros(1,length(x));
    y(find(cell_labels==3)) = 1;
        
    figure; plot(x,y,'ko')
    ylim([-0.2 1.2])
    set(gca,'ytick',[0,1])
    xlabel('Marker gene expression');
    ylabel('Beta cell or not');
end

if (part_2)
    b = glmfit(x', y', 'binomial', 'link', 'logit');
    yfit = glmval(b, x, 'logit');

    figure; plot(x,y,'ko')
    ylim([-0.2 1.2]);     set(gca,'ytick',[0,1]);    xlabel('Marker gene expression');     ylabel('Beta cell or not');
    hold on; plot(x, yfit, 'ro');
    
%     1 / (1+ exp(-1*(x(1)*b(2)+b(1))))
end

if (part_3)

    rng(1)
    r = randperm(500);
    
    %training set
    x_training = x(r(1:300));
    y_training = y(r(1:300));

    %test set
    x_test = x(r(301:500));
    y_test = y(r(301:500));
    
    b_training = glmfit(x_training', y_training', 'binomial', 'link', 'logit');

    figure; plot(x_test,y_test,'ko')
    y_test_fit = glmval(b_training, x_test, 'logit');

    ylim([-0.2 1.2]);     set(gca,'ytick',[0,1]);    xlabel('Marker gene expression');     ylabel('Beta cell or not');
    hold on; plot(x_test, y_test_fit, 'ro');

    threshold = 0.5;
    y_test_predict = y_test_fit;
    y_test_predict(find(y_test_predict<threshold))  = 0;
    y_test_predict(find(y_test_predict>=threshold)) = 1;
    crosstab(y_test', y_test_predict)
end

