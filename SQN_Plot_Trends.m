%% read data from Trend file

[filename,path] = uigetfile('*.tdms','Select file to load');
fname1=fullfile(path,filename);
[time,nRMS,nPK,nCF,vPK,dPK]=VSreadtrend(fname1,'Native RMS','Native 0-pk','Native pk/RMS','Sgl-Int 0-pk','Dbl-Int 0-pk');

%% read plant parameter data from TXT file

[filename,path] = uigetfile('*.txt','Select file to load');
fname2=fullfile(path,filename);
ss=fileread(fname2);
dd=textscan(ss,'%s %s %f %f %f');
temp=strcat(dd{1,1},{' '},dd{1,2});
pptime=datetime(temp);
ppdata=zeros(length(pptime),3);
for i=1:3
    ppdata(:,i)=dd{1,2+i};
end

%% read dataset list from TXT file

if ~exist('fname3','var')
    [filename,path] = uigetfile('*.txt','Select file to load');
    fname3=fullfile(path,filename);
end
ss=fileread(fname3);
dd=textscan(ss,'%s %s',1,'Delimiter','\t');
fpath=char(dd{1,2});
dd=textscan(ss,'%s %s %s','headerLines',3,'Delimiter','\t');
for i=1:length(dd{1,1})
    wfiles(i)=fullfile(fpath,dd{1,1}(i,1));
    tstart=datetime(dd{1,2}(i,1),'InputFormat','MM/dd/yy HH:mm:ss');
    tend=datetime(dd{1,3}(i,1),'InputFormat','MM/dd/yy HH:mm:ss');
    dstimes{i}=[tstart tend];
end

%% plot various trends

ts=datetime(2020,05,06,20,00,00);
te=datetime(2020,05,09,00,00,00);

%acceleration RMS
ymax=80;
h=figure();
t=tiledlayout(2,2,'TileSpacing','compact');
for i=1:4
    nexttile;
    hold on;
    for j=1:length(dstimes)
        fx(j)=fill([dstimes{j}(1) dstimes{j}(1) dstimes{j}(2) dstimes{j}(2)],[0 ymax ymax 0],'r','LineStyle','none');
        tmid=dstimes{j}(1)+(dstimes{j}(2)-dstimes{j}(1))/2;
        text(tmid,ymax,char(64+j),'HorizontalAlignment','center','VerticalAlignment','top','Color','r');
    end
    for j=1:3
        plot(time,nRMS(:,(i-1)*3+j));
    end
    hold off;
    ylim([0 ymax]);
    ylabel('Acceleration (g, RMS)');
    yyaxis right;
    plot(pptime,ppdata(:,1),'--k')
    ylim([0 100]);
    ylabel('Reactor Power (%)');
    ytickformat('percentage');
    xlim([ts te]);
    title(strcat('Sensor-',num2str(i)));
    for j=1:length(dstimes)
        set(get(get(fx(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    legend('X-Dir','Y-Dir','Z-Dir','%RXP','Location','northeast');
    ax=gca;
    ax.YAxis(2).Color='k';
end
title(t,'Acceleration Trends (RMS)');

%acceleration peak
ymax=200;
h=figure();
t=tiledlayout(2,2,'TileSpacing','compact');
for i=1:4
    nexttile;
    hold on;
    for j=1:length(dstimes)
        fx(j)=fill([dstimes{j}(1) dstimes{j}(1) dstimes{j}(2) dstimes{j}(2)],[0 ymax ymax 0],'r','LineStyle','none');
        tmid=dstimes{j}(1)+(dstimes{j}(2)-dstimes{j}(1))/2;
        text(tmid,ymax,char(64+j),'HorizontalAlignment','center','VerticalAlignment','top','Color','r');
    end
    for j=1:3
        plot(time,nPK(:,(i-1)*3+j));
    end
    hold off;
    ylim([0 ymax]);
    ylabel('Acceleration (g, 0-pk)');
    yyaxis right;
    plot(pptime,ppdata(:,1),'--k')
    ylim([0 100]);
    ylabel('Reactor Power (%)');
    ytickformat('percentage');
    xlim([ts te]);
    title(strcat('Sensor-',num2str(i)));
    for j=1:length(dstimes)
        set(get(get(fx(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    legend('X-Dir','Y-Dir','Z-Dir','%RXP','Location','northeast');
    ax=gca;
    ax.YAxis(2).Color='k';
end
title(t,'Acceleration Trends (Peak)');

%acceleration crest factor
ymax=20;
h=figure();
t=tiledlayout(2,2,'TileSpacing','compact');
for i=1:4
    nexttile;
    hold on;
    for j=1:length(dstimes)
        fx(j)=fill([dstimes{j}(1) dstimes{j}(1) dstimes{j}(2) dstimes{j}(2)],[0 ymax ymax 0],'r','LineStyle','none');
        tmid=dstimes{j}(1)+(dstimes{j}(2)-dstimes{j}(1))/2;
        text(tmid,ymax,char(64+j),'HorizontalAlignment','center','VerticalAlignment','top','Color','r');
    end
    for j=1:3
        plot(time,nCF(:,(i-1)*3+j));
    end
    hold off;
    ylim([0 ymax]);
    ylabel('Acceleration (g, pk/RMS)');
    yyaxis right;
    plot(pptime,ppdata(:,1),'--k')
    ylim([0 100]);
    ylabel('Reactor Power (%)');
    ytickformat('percentage');
    xlim([ts te]);
    title(strcat('Sensor-',num2str(i)));
    for j=1:length(dstimes)
        set(get(get(fx(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    legend('X-Dir','Y-Dir','Z-Dir','%RXP','Location','northeast');
    ax=gca;
    ax.YAxis(2).Color='k';
end
title(t,'Acceleration Trends (Crest Factor)');

%velocity peak
ymax=10;
h=figure();
t=tiledlayout(2,2,'TileSpacing','compact');
for i=1:4
    nexttile;
    hold on;
    for j=1:length(dstimes)
        fx(j)=fill([dstimes{j}(1) dstimes{j}(1) dstimes{j}(2) dstimes{j}(2)],[0 ymax ymax 0],'r','LineStyle','none');
        tmid=dstimes{j}(1)+(dstimes{j}(2)-dstimes{j}(1))/2;
        text(tmid,ymax,char(64+j),'HorizontalAlignment','center','VerticalAlignment','top','Color','r');
    end
    for j=1:3
        plot(time,vPK(:,(i-1)*3+j));
    end
    hold off;
    ylim([0 ymax]);
    ylabel('Velocity (ips, 0-pk)');
    yyaxis right;
    plot(pptime,ppdata(:,1),'--k')
    ylim([0 100]);
    ylabel('Reactor Power (%)');
    ytickformat('percentage');
    xlim([ts te]);
    title(strcat('Sensor-',num2str(i)));
    for j=1:length(dstimes)
        set(get(get(fx(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    legend('X-Dir','Y-Dir','Z-Dir','%RXP','Location','northeast');
    ax=gca;
    ax.YAxis(2).Color='k';
end
title(t,'Velocity Trends (Peak)');


%% call batchTDMSrun function to process each selected dataset

[dsets]=batchTDMSrun('FileList',wfiles,'FileSubsets',dstimes,'StoreWaveforms',true,'CalculateSpectra',true);

%% TESTING/VALIDATION CODE - build 'finfo' structure containing parameters associated with each file

for i=1:length(wfiles)
    
    %specify folder and filename, open using "no data" option for max speed
    [fpath,fname,fext]=fileparts(wfiles{i});
    finfo(i).folder=fpath;
    finfo(i).name=strcat(fname,fext);
    [~,metastruct]=TDMS_readTDMSFile(fullfile(finfo(i).folder,finfo(i).name),'UTC_DIFF',0,'GET_DATA_OPTION','getnone');
    
    %determine/set UTC offset parameter
    idx=find(strcmp(metastruct.rawDataInfo(1).propNames,'UTC Offset'));
    if ~isempty(idx)
        finfo(i).utc=metastruct.rawDataInfo(1).propValues{idx};
    else
        finfo(i).utc=0;
    end
    
    %determine the sampling rate and number of points
    idx=find(strcmp(metastruct.rawDataInfo(3).propNames,'wf_increment'));
    finfo(i).fs=double(1/metastruct.rawDataInfo(3).propValues{idx});
    finfo(i).numpts=metastruct.numberDataPoints(3);
    
    %determine the start and end times for the selected file
    idx=find(strcmp(metastruct.rawDataInfo(3).propNames,'wf_start_time'));
    sstart=metastruct.rawDataInfo(3).propValues{idx};
    finfo(i).tstart=datetime(datetime(sstart,'InputFormat','dd-MMM-yyyy HH:mm:ss:SSS')+hours(dt),'Format','dd-MMM-yyyy HH:mm:ss:SSS');
    finfo(i).tend=finfo(i).tstart+seconds((finfo(i).numpts-1)/finfo(i).fs);
    
    %if data subsets are to be selected, determine if they are in-range
    if ~isempty(fsubsets);
        fbounds=fsubsets{i};
        finfo(i).inrange=0;
        if length(fbounds)==2
            sval=fbounds(1);
            eval=fbounds(2);
            if isdatetime(sval) %evaluate start time first
                sflag=isbetween(sval,finfo(i).tstart,finfo(i).tend);
                is=round(seconds(sval-finfo(i).tstart)*finfo(i).fs)+1;
            else
                sflag=sval>1&&sval<finfo(i).numpts;
                is=sval;
            end
            if isdatetime(sval) %evaluate end time next
                eflag=isbetween(eval,finfo(i).tstart,finfo(i).tend);
                ie=round(seconds(eval-finfo(i).tstart)*finfo(i).fs)+1;
            else
                eflag=eval>1&&eval<finfo(i).numpts;
                ie=eval;
            end
            if sflag && eflag
                finfo(i).inrange=1;
                finfo(i).istart=is;
                finfo(i).iend=ie;
            end
        end
    end
    
end

%% TESTING/VALIDATION - load a subset of a specified TDMS file (verifying functionality)

stemp=datetime(2020,03,27,14,09,00);
etemp=datetime(2020,03,27,14,10,00);

ftemp='C:\SI Projects\1901462 (BFN Valves)\Data\Startup\01a\Auto_2020_03_27_120822_Wave.tdms';
[~,metastruct]=TDMS_readTDMSFile(ftemp,'UTC_DIFF',0);

idx=find(strcmp(metastruct.rawDataInfo(1).propNames,'UTC Offset'));
temp=metastruct.rawDataInfo(1).propValues(idx);
dt=temp{1};

ist = find(strcmp(metastruct.rawDataInfo(3).propNames,'wf_start_time'));
iinc = find(strcmp(metastruct.rawDataInfo(3).propNames,'wf_increment'));
startstr = metastruct.rawDataInfo(3).propValues{ist};
numPts = metastruct.numberDataPoints(3);
Fs=double(1/metastruct.rawDataInfo(3).propValues{iinc});
tstart = datetime(datetime(startstr,'InputFormat','dd-MMM-yyyy HH:mm:ss:SSS')+hours(dt),'Format','dd-MMM-yyyy HH:mm:ss:SSS');
tend = tstart+seconds((numPts-1)/Fs);

istart = round(seconds(stemp-tstart)*Fs)+1;
iend = round(seconds(etemp-tstart)*Fs);

data = TDMS_readTDMSFile(ftemp,'SUBSET_GET',[istart iend],'SUBSET_IS_LENGTH',false);
