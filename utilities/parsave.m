% save data (necessary for use within parallel process)
function parsave(data,param,fname,ll,InitParam,valLL)
if nargin < 3
    save(param,'-v7.3','data');
elseif nargin < 4
    save(fname,'-v7.3','data','param');
elseif nargin <5
    save(fname,'-v7.3','data','param','ll');
elseif nargin < 6
    save(fname,'-v7.3','data','param','ll','InitParam');
else
    save(fname,'-v7.3','data','param','ll','InitParam','valLL');
end
end