%% Function BINDATA
%  Bin the data with binwidth specified in Param
%
% Input:
%       Param - Structure containing commonly used parameters
%       data  - Matrix containing data to be binned along 2nd dimension
%                   dim(NxT)
%
% Output:
%       binnedData - Binned data
%                       dim(NxT/BIN)
%
%% AUTHOR    : Koh Tze Hui
%% DEVELOPED : MATLAB (R2017a)
%% FILENAME  : bindata.m


function binnedData = bindata(Param,data,avg)
if nargin<3
    avg = Param.BIN; %average over bin
end
[nData,tData] = size(data);
y = reshape(data,nData,Param.BIN,tData/Param.BIN);
binnedData = squeeze(sum(y,2))/avg;
end