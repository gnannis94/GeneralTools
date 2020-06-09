function [time,varargout] = VSreadtrend(filename,varargin)
%testting testing testingS
%look at me i'mbreaking everything
%hey break this break that yea yea yea
%ensure that the TDMS_reader functions are added to the search path
spaths=regexp(path,pathsep,'split');
idx=find(contains(spaths,'TDMS_reader'));
if isempty(idx)
    [startpath,~,~]=fileparts(mfilename('fullpath'));
    tgtpath=fullfile(startpath,'TDMS_reader');
    if exist(tgtpath,'dir')
        addpath(genpath(tgtpath));
    else
        error('Could not locate the "TDMS_reader" function necessary to process the file.');
    end
end

filepath='C:\SI Software\Matlab\Miscellaneous\TDMS Read\TDMS_reader\';
addpath(genpath(filepath));

%read partial TDMS file (meta infomration only) to determine whether UTC
%offset parameter exists, and take action accordingly
[~,meta]=TDMS_readTDMSFile(filename,'GET_DATA_OPTION','getnone');
idx=find(strcmp(meta.rawDataInfo(1).propNames,'UTC Offset'));
if isempty(idx)
    local=hours(tzoffset(datetime('now','TimeZone','local')));
    prompt=['Enter Desired UTC Offset (default=0)' sprintf('\n') '[Note: your local time is UTC ' num2str(local) ']'];
    answer=inputdlg(prompt,'Specify UTC Offset',[1 45],{'0'});
    dt=str2num(answer{1});
else
    temp=meta.rawDataInfo(1).propValues(idx);
    dt=temp{1};
end

% read full TDMS file from disk (time-intensive, could be modified in the future
% to only read user-specified values)
[output,meta]=TDMS_readTDMSFile(filename,'UTC_DIFF',0);

% read the datetime information for the TDMS file by finding the index of the appropriate group
% adjust the time zone of the data based on the read/specified UTC offset
itime=find(strcmp(meta.chanNames,'Time'));
numpts=meta.numberDataPoints(itime);
time=datetime(output.data{1,itime},'ConvertFrom','datenum')+hours(dt);
time=dateshift(time,'start','second','nearest').';

% determine the channel indices for each specified datagroup (varargin), initialize,
% and write a 2D array containing those values to varargout
for i=1:nargin-1
    igroup=find(strcmp(output.groupNames,varargin{i}));
    ichans=output.chanIndices{:,igroup};
    icmin=min(ichans);
    icmax=max(ichans);
    numchs=icmax-icmin+1;
    data=zeros(numpts,numchs);
    for j=1:numchs
        data(:,j)=output.data{1,icmin+j-1};
    end
    varargout{i}=data;
end

end