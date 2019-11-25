% dfc-lrr source code:
% https://www.cs.cmu.edu/~atalwalk/dfc_lrr_matlab.tar.gz
% clc;
clear;
close all;

%% general setting
%#ok<*NBRAK>
if ~exist('./Output','dir')
   mkdir('./Output');
end
addpath ('../common_usage/');
addpath('../data/UCSD/')
addpath('../data/COIL20/')
addpath('../data/ExtYaleB/')
addpath('../data/MNIST/')

%% load data
% V5_32_32_COIL20
% V20_72_32_32_COIL20
% V20_36_32_32_COIL20

% V21-38_48_42_ExtYaleB
% V1-20_48_42_ExtYaleB

% V3_90_135_UCSD
% V10_10000_20_20_MNIST_Test

clear
par.dataSet = 'V21-38_48_42_ExtYaleB';
par.lambda = [0.1]; %#ok<*NBRAK>
par.post = 1;
par.normalizeMethod = 'max';
par.part = [1];
par.run = 1;
main_dfc_lrr(par)