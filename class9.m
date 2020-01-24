load thyroid_working.mat A features pathology
% figure; hist(A(:,2))

response = strmatch('bestresp',features);
age      = strmatch('agest',   features);

% figure; plot(A(:,age),A(:,response),'o');

older_indicator = zeros(1,size(A,1));
older_indicator(find(A(:,age)> median(A(:,2))) ) = 1; 


%lets do a t-test
%response of the younger versus the response of the older

[H, P] = ttest2(A(find(older_indicator==0),response), A(find(older_indicator==1),response));

%let's make a boxplot
M = nan(size(A,1),2);
M(find(older_indicator==0),1) = A(find(older_indicator==0),response);
M(find(older_indicator==1),2) = A(find(older_indicator==1),response);
% figure; boxplot(M)

%lets look at a discrete variable

hand_foot = strmatch('hfsyn',   features);
hand_foot_indicator = zeros(1,size(A,1));
hand_foot_indicator(find(A(:,hand_foot)> 1) ) = 1; 
[TABLE,CHI2,P,LABELS] = crosstab(older_indicator, hand_foot_indicator);

%can we use a non-parametric test? (if we can't assume a normal
%distribution)

[P, H] = ranksum(A(find(older_indicator==0),response), A(find(older_indicator==1),response))


% If time allows:
% corr([A(:,age),A(:,response)])


% figure; hist(A(:,response),6)
% xlabel('Tumor reduction'); ylabel('Frequency');

