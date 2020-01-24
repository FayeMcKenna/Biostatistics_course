%% clean up
clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same

[A, names] = xlsread('C_elegans_dataset.xlsx'); %read in dataset

A_mat = A(:,3:end); % slice A to keep only numeric values

gene_names = names(2:end, 1); % save gene names without header

A = xlsread("C_elegans_dataset.xlsx", 1, "C2:AZ20518"); % the second input "1" is to specify which sheet of the excel file to import, can also leave as empty, i.e. ""

gene_names = xlsread("C_elegans_dataset.xlsx", 1, "A2:A20518");

figure; imagesc(A); %image as matrix
colorbar
