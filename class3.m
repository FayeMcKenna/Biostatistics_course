% PCA- reduces dimensionality 

%% clean up
clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same
% one dimension shows the most variation in one way and the other dimension
% the least variation, look at broad patterns and use those
%load pancreas data set
try
    A;
catch
    load pancreas_data.mat
end    


A_for_pca=log10(1+A);
[A_COEFF, A_SCORE,latent,tsquared,explained,mu]=pca(A_for_pca');



%% class 3

X= A_for pca