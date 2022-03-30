%% Function SPLITTRIAL
%  Takes in a concatenated matrix and splits it back into trials
%
%  Input:
%       Data      - Structure containing data to be split into trials
%                       
%       N         - Scalar number of trials in matrix
%       fieldName_1 - String containing name of field to be split
%       fieldName_2 - String containing name of field to be split
%
%  Output:
%       Z         - Structure containing data split by trials with
%                   fieldnames fieldName_1 and fieldnName_2
%                       dim(1xN_TRIAL)
%
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2018a)
%% FILENAME  : splittrial.m

function Z = splittrial(Data,N,fieldName_1,fieldName_2)
if nargin < 4
    fieldName_2 = nan;
end

if isstruct(Data)
z = [Data(:).(fieldName_1)];
if nargin >3
x = [Data(:).(fieldName_2)];
end
else
    z = Data;
end
triallen = size(z,2)/N;
Z = struct;
for i = 1:N
    Z(i).(fieldName_1) = z(:,((i-1)*triallen+1):triallen*i);
    if nargin > 3
    Z(i).(fieldName_2) = x(:,((i-1)*triallen+1):triallen*i);
    end
end
end