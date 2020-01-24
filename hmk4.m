%% Homwork 3

clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same

%% Question 1 

nump=23; 
samples=10000; 
birthd=ceil(365*rand(nump,samples)); 
count=0; 
for j=1:samples 
if numel(birthd(:,j))-numel(unique(birthd(:,j))) >0 
count=count+1; 
end 
end 
prob= 1-count/samples
prob2 = (1 - prob )*100
% prob3 = 100 - prob2


%% Question 2
% Get the values for the domain, x.
x = 0:10;
z = 0:100;
% Get the values of the probability mass function. 
a = binopdf(x, 10, 0.1);
b = binopdf(x, 10, 0.6);
c = binopdf(z, 100, 0.1);
d = binopdf(z, 100, 0.6);

% Do the plots. 
subplot(1,4,1),bar(x,a,1,'b') 
title(' n = 10, p = 0.1') 
xlabel('X'),ylabel('f(X)')
axis square
subplot(1,4,2),bar(x,b,1,'g')
title(' n = 10, p = 0.6') 
xlabel('X'),ylabel('f(X)')
axis square
subplot(1,4,3),bar(z,c,1,'r')
title(' n = 100, p = 0.1') 
xlabel('X'),ylabel('f(X)')
axis square
subplot(1,4,4),bar(z,d,1,'y')
title(' n = 100, p = 0.6')
xlabel('X'),ylabel('f(X)')
axis square

%% Question 3: 

% for ii = 0:99
%    matrix{ii+1} = randi([0 1],1000,1);
% end


% create 100 random vectors with only zeros and ones each having a length
% 1000

r = randi([0 1],1000,100);
a_diff = randi([0 1],1000,100);

% identify longest stretches of each of the vectors
for j = 1:100
    a = r(:,j);
     for c = 1:100
        a_diff = diff(a);
        break_idx = find([1; a_diff; 1] ~= 0);
        stretches = diff(break_idx);
        longest_stretch = max(stretches);
     end
     longstat(j,:)=[longest_stretch,0];
end

%plot a histogram of lengths

A = longstat(:,1);
scaled_A= log10(1+A);index=scaled_A~=0;

m=figure; colorbar;hist(scaled_A, 10);
xlabel('vectors');
ylabel('longest stretch');
title('Distribution of Longest Stretch in 100 Vectors');
legend('mean: 10.4300 median: 10');% legend

% mean and median values of plot

[x y]=imhist(m);
count=(sum(x(:)))/2;


mean_value= mean(A)
median_value= median(A)

% create 100 random vectors with only zeros and ones each having a length
% 100

r = randi([0 1],100,100);
a_diff = randi([0 1],1000,100);

% identify longest stretches of each of the vectors
for j = 1:100
    a = r(:,j);
     for c = 1:100
        a_diff = diff(a);
        break_idx = find([1; a_diff; 1] ~= 0);
        stretches = diff(break_idx);
        longest_stretch = max(stretches);
     end
     longstat(j,:)=[longest_stretch,0];
end

%plot a histogram of lengths

A = longstat(:,1);
scaled_A= log10(1+A);index=scaled_A~=0;

figure; colorbar;hist(scaled_A, 10);
xlabel('vectors');
ylabel('longest stretch');
title('Distribution of Longest Stretch in 100 Vectors');
legend('mean: 6.8 median: 6.5');% legend

% mean and median values of plot
mean_value= mean(A)
median_value= median(A)
%% Question 4

load('thyroid_working.mat')

%look at dose
dose = A(:,14);

% take out NAN
k = find(isnan(dose))'; 
dose(k) = 0;
dose(isnan(dose)) = 0;

% determine cut off 650
mean_dose= mean(dose); % 554.5455
median_dose= median(dose); %600

% make groups based on cut off
x = dose;
nCol=1; 
x(x(:,nCol)==0)=nan;% say, for the particular column
x(x(:,nCol)<650)=1;
x(x(:,nCol)>650)=2;


% match with response
bestresp = A(:,19);


% add together variables and sort based on groups
test= [x,bestresp];

B = sortrows(test);

low = B(1:22,2);
high = B(23:44,2);

lowhigh = [low,high];

% two-samples t test
x = lowhigh(:,1); % groups
y = lowhigh(:,2); % score


%testing that the mean of high group is higher than low group
[h,p] = ttest2(x,y,'Tail','left','Alpha',0.05,'Vartype','equal');

% p = .96
% h=0 does not reject the null hypothesis



