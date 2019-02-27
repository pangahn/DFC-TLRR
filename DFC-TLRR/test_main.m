clear
close all

%% addpath
current_path = fileparts(mfilename('fullpath'));
addpath(genpath(current_path));
addpath(genpath('../data'));
addpath ('../common_usage/')
%addpath ('../common_usage/tensor_toolbox/')    
    
%% general setting
%#ok<*NBRAK>
if ~exist('./Output','dir')
   mkdir('Output');
end

%% optional data sets
% T20_72_32_32_COIL20
% T5_32_32_COIL20
% T100_72_32_32_COIL100
% T100_36_32_32_COIL100

% T10_48_42_ExtYaleB
% T38_48_42_ExtYaleB
% T1-10_48_42_ExtYaleB
% T11-20_48_42_ExtYaleB
% T21-30_48_42_ExtYaleB
% V31-38_48_42_ExtYaleB

% T3_90_135_UCSD

%% serialization
clear
parms.display = true;
parms.dataSet = 'T20_72_32_32_COIL20';
parms.normalizeMethod = 'l2';
parms.method = 'nuclear';
parms.lambda = [4];
parms.partition = [5];
parms.s = [0.99];
parms.run = 1;
main(parms)
