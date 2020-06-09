function [dataStruct] = batchTDMSrun(varargin);
% Read all TDMS files in user-specified directory and process with a
% consistent set of parameters.

%parse incoming function inputs to determine whether
p=inputParser;
addParameter(p,'FileList','',@iscell);
addParameter(p,'FileSubsets','',@ismatrix);
addParameter(p,'MainPath','',@ischar);
addParameter(p,'SearchSubfolders',false);
addParameter(p,'VirtualCombs',@iscell);
addParameter(p,'StoreWaveforms',false);
addParameter(p,'CalculateSpectra',false);
parse(p,varargin{:});

%ensure that the TDMS_reader functions are added to the search path
spaths=regexp(path,pathsep,'split');
idx=find(contains(spaths,'TDMS_reader'));
if isempty(idx)
    [startpath,~,~]=fileparts(mfilename,'fullpath');
    tgtpath=fullfile(startpath,'TDMS_reader');
    if exist(tgtpath,'dir')
        addpath(genpath(tgtpath));
    else
        error('Could not locate the "TDMS_reader" function necessary to process the file.');
    end
end

% file selection switch, based on the variables that were passed
if ~isempty(p.Results.FileList)
    for i=1:length(p.Results.FileList)
        [fpath,fname,fext]=fileparts(p.Results.FileList{i});
        dirFiles(i).folder=fpath;
        dirFiles(i).name=strcat(fname,fext);
    end
    [dirFiles(1).folder '\' dirFiles(1).name]
else
    % Provides user dialog to select TDMS folder of interest (if not passed)
    if ~isempty(p.Results.MainPath)
        fpath=p.Results.MainPath;
    else
        fpath = uigetdir(pwd,'Select folder with data files:');
    end
    
    % Finds all _Wave TDMS files in indicate path (or subfolders, if that
    % option was selected)
    if p.Results.SearchSubfolders
        dirFiles = dir([fpath '\**\*_Wave.tdms']);
    else
        dirFiles = dir([fpath '\*_Wave.tdms']);
    end
end

[~,metastruct] = TDMS_readTDMSFile([dirFiles(1).folder '\' dirFiles(1).name],'GET_DATA_OPTION','getnone');

incInd = find(strcmp(metastruct.rawDataInfo(3).propNames,'wf_increment'));    
fs = double(1/metastruct.rawDataInfo(3).propValues{incInd});

[h ,param] = ProcessingParameters(fs);

% Bandpass Filter parameters
bandflag = param.bandtog;
n = param.nB/2;
fLow = param.flow;
fHigh = param.fhigh;
ripple = param.ripple;

% Notch Filter parameters
notchflag = param.notchtog;
nN = param.nN;
filtflag = param.notchflag;
notchfreq = param.notchfreq;
notchwidth = param.notchwidth;

% Integration parameters
intflag = param.inttog;
highflag = param.hightog;
nH = param.nH;
fc = param.fc;

% waveform trucation
wavetrunc = param.wavetrunc;

% frequency spectra computation (placeholder - need form fields)
fcut=0;
fmaxN=3000;
fmaxV=300;
fmaxD=50;
dF=0.1;
ovlp=90;

w = waitbar(0,'Processing In Progress...','Name','Progress');
pause(2)

for j = 1:length(dirFiles)

waitbar(((j-1))/length(dirFiles),w,['Loading Data File ' num2str(j) '...'])     
pause(2)

enum=j;
filename = [dirFiles(j).folder '\' dirFiles(j).name];
dataStruct(enum).filepath = dirFiles(j).folder;
dataStruct(enum).filename = dirFiles(j).name;

% First call to TDMS_readTDMSFile does not upload data. Meta structure is
% obtained to retrieve information about the data in the TDMS file
[~,metastruct] = TDMS_readTDMSFile(filename,'GET_DATA_OPTION','getnone');

numPts = metastruct.numberDataPoints(3);
incInd = find(strcmp(metastruct.rawDataInfo(3).propNames,'wf_increment'));    
fs = double(1/metastruct.rawDataInfo(3).propValues{incInd});
fnyq = fs/2;

trunc = floor(wavetrunc*fs);

chanN = metastruct.chanNames(3:length(metastruct.chanNames));

data = TDMS_readTDMSFile(filename); 

dataStruct(enum).numPts = numPts;

numImport = length(data.data{3});

dataStruct(enum).numImport = numImport;

%read datarray from imported data struct
for i=3:length(data.data)
    datarray(:,i-2) = double(data.data{1,i}');
end

%create virtual channels if requested/prompted
if ~isempty(p.Results.VirtualCombs)
    vchans=p.Results.VirtualCombs;
    for i=1:length(vchans)
        temp=zeros(length(datarray),1);
        for k=1:length(vchans{i})
            temp=temp+datarray(:,vchans{i}(k)); %sum specified channels
        end
        datarray=[datarray temp/length(vchans{i})]; % append average to fdata
    end
end

waitbar(((j-1)+1/3)/length(dirFiles),w,['Data File ' num2str(j) ' Loaded'])     
pause(1)

waitbar(((j-1)+1/4)/length(dirFiles),w,['Filtering Data File ' num2str(j) '...'])     
pause(1)

%bandpass filtering
for i=1:size(datarray,2)
    if bandflag == 1
        [numband,denband] = cheby1(n,ripple,[fLow fHigh]/fnyq,'bandpass');
        datarrayF(:,i) = filtfilt(numband,denband,datarray(:,i));
        
        dataStruct(enum).n = n;
        dataStruct(enum).fLow = fLow;
        dataStruct(enum).fHigh = fHigh;
    else
        datarrayF(:,i) = datarray(:,i);
    end
end

%notch filtering
if notchflag == 1
    nind = find(filtflag == 1);
    if numel(nind) > 0
        for k=1:length(nind)
            f1 = notchfreq(nind(k))-notchwidth(nind(k))/2;
            f2 = notchfreq(nind(k))+notchwidth(nind(k))/2;
            
            [numnotch,dennotch] = cheby1(nN,ripple,[f1 f2]/fnyq,'stop');
            
            datarrayF = filtfilt(numnotch,dennotch,datarrayF);
            
            dataStruct(enum).notchfreq(k) = notchfreq(nind(k));
            dataStruct(enum).notchwidth(k) = notchwidth(nind(k));
        end
        
        dataStruct(enum).nN = nN;
    end
end

%subtract mean
datarrayF = datarrayF-mean(datarrayF);
sizeDat = size(datarrayF);

if intflag == 1
    
    %digital integration
    waitbar(((j-1)+2/4)/length(dirFiles),w,['Integrating Data File ' num2str(j) '...'])     
    pause(1)
    [numhigh,denhigh] = cheby1(nH,ripple,fc/fnyq,'high');

    Vint = 386.4*cumsum(datarrayF)/fs;
    if highflag == 1
        dataStruct(enum).nH = nH;
        dataStruct(enum).fc = fc;
        velo = filtfilt(numhigh,denhigh,Vint);
    else
        velo = Vint;
    end
    velocut = velo(trunc+1:numImport-trunc,:)-mean(velo(trunc+1:numImport-trunc,:));
    
    Dint = 1000*cumsum(velo)/fs;
    if highflag == 1
        disp = filtfilt(numhigh,denhigh,Dint);
    else
        disp = Dint;
    end
    dispcut = disp(trunc+1:numImport-trunc,:)-mean(disp(trunc+1:numImport-trunc,:));
        
else
    dataStruct(enum).velocity = zeros(sizeDat(1),sizeDat(2));
    dataStruct(enum).displacement = zeros(sizeDat(1),sizeDat(2));
end

%truncate processed nativewaveform
datarrayF(trunc+1:numImport-trunc,:);

%calculate frequency spectra (if enabled)
if p.Results.CalculateSpectra
    
    waitbar(((j-1)+3/4)/length(dirFiles),w,['Calculating PSD/Spectra for Data File ' num2str(j) '...'])     
    pause(1)
    
    BSavg=2*floor(fs/dF/2);
    dFavg=fs/BSavg;
    noverlap=round(BSavg*ovlp/100);
    
    icut=floor(fcut/dFavg)+1;
    imaxN=floor(fmaxN/dFavg)+1;
    imaxV=floor(fmaxV/dFavg)+1;
    imaxD=floor(fmaxD/dFavg)+1;
    
    [psdN,freqN]=pwelch(datarrayF,hann(BSavg),noverlap,BSavg,fs);
    psdN(1:icut,:)=0;
    psdN=psdN(1:imaxN,:);
    dataStruct(enum).psdN=psdN;
    dataStruct(enum).spectN=sqrt(1.5*psdN*dFavg);
    dataStruct(enum).freqN=freqN(1:imaxN);
    
    [psdV,freqV]=pwelch(velocut,hann(BSavg),noverlap,BSavg,fs);
    psdV(1:icut,:)=0;
    psdV=psdV(1:imaxV,:);
    dataStruct(enum).psdV=psdV;
    dataStruct(enum).spectV=sqrt(1.5*psdV*dFavg);
    dataStruct(enum).freqV=freqV(1:imaxV);
    
    [psdD,freqD]=pwelch(dispcut,hann(BSavg),noverlap,BSavg,fs);
    psdD(1:icut,:)=0;
    psdD=psdD(1:imaxD,:);
    dataStruct(enum).psdD=psdD;
    dataStruct(enum).spectD=sqrt(1.5*psdD*dFavg);
    dataStruct(enum).freqD=freqD(1:imaxD);
    
    dataStruct(enum).dFavg=dFavg;
    dataStruct(enum).ovlp=ovlp;
    
end

% Calculates the maximum (peak) values
dataStruct(enum).Nmax = max(abs(datarrayF));
dataStruct(enum).Vmax = max(abs(velocut));
dataStruct(enum).Dmax = max(abs(dispcut));

% Calculating RMS and 95% confidence amplitudes
for h = 1:sizeDat(2)
    dataStruct(enum).Nrms(h) = norm(datarrayF(:,h))/sqrt(sizeDat(1));
    dataStruct(enum).Nconf95(h) = 1.96*dataStruct(enum).Nrms(h);
    dataStruct(enum).Vrms(h) = norm(velocut(:,h))/sqrt(sizeDat(1));
    dataStruct(enum).Vconf95(h) = mean(velocut(:,h))+1.96*dataStruct(enum).Vrms(h);
    dataStruct(enum).Drms(h) = norm(dispcut(:,h))/sqrt(sizeDat(1));
    dataStruct(enum).Dconf95(h) = mean(dispcut(:,h))+1.96*dataStruct(enum).Drms(h);
end

% Storing waveforms in data structure (if enabled)
if p.Results.StoreWaveforms
    dataStruct(enum).dataRnat = datarray;
    dataStruct(enum).timeR = (0:length(dataStruct(enum).dataNr)-1)/fs;
    dataStruct(enum).dataPnat = datarrayF;
    dataStruct(enum).timeP = (0:length(dataStruct(enum).dataNp)-1)/fs;
    dataStruct(enum).dataPvel = velocut;
    dataStruct(enum).dataPdisp = dispcut;
end

dataStruct(enum).fs = fs;
dataStruct(enum).ripple = ripple;
dataStruct(enum).wavetrunc = wavetrunc;

clear datarray datarrayF fs

clear sizeDat chani dashind Dint disp dispcut i Vint velo velocut filename...
      enum numband denband numhigh denhigh data fnyq 

waitbar(((j-1)+1)/length(dirFiles),w,['Data File ' num2str(j) ' Processed'])     
pause(2)
  
end
 
waitbar(1,w,'Complete!')     
pause(1)
close(w)

clear bandflag chanN dennotch dirFiles f1 f2 fc fHigh fLow filtflag h...
      highflag idx incInd intflag j k metastruct n nH nind nN notchflag...
      notchfreq notchwidth numImport numnotch numPts path ripple spaths...
      trunc w wavetrunc
  
end