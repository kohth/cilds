%% Function CONCATENATE
%  Concatenates structure field into dim(M, NxT)
%
%  Input:
%       X    - 1xN struct 
%                   1. trialId - scalar idenfication of trial
%                   2. 'fieldName' - Field to concatenate
%                          dim(M,T)
%
%       fieldName - Identifier of field to be concatentated
%
%  Output:
%       z - Matrix of concatenated data
%               dim(M, NxT)
%
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2017a)
%% FILENAME  : concatenate.m

function z = concatenate(X,fieldName)
N = size(X,2);
% [M,T] = size(X(1).(fieldName));
z = X(1).(fieldName);
for i = 2:N
    z = horzcat(z,X(i).(fieldName));
end
end