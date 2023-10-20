% convert uManTime to stimTime

imFile='C:\Users\jshha\Dropbox\ColonRamosLab\CaImaging_uMan\181009_AFDeat4_DCR7521\181009_AFDeat4_DCR7521_MMStack_Pos0.ome.tif';
% insert guiSelect vid file and strip rest of the locations...

imDir=fileparts(imFile);
baseDir=fileparts(imDir);
stimDir=fullfile(baseDir,'stimFiles');


[vTimeSec,elapsedT,iDay] = getuManTime_Yale('imFull',imFile)

d=strcat(iDay(6:7),iDay(9:10),iDay(3:4));
hT=vTimeSec(1)/(60*60);
h=floor(hT);
mT=(hT-h)*60;
m=floor(mT);
sT=(mT-m)*60;
s=floor(sT);

dStr=strcat(d,num2str(h),num2str(m),num2str(s));
vStr=str2num(dStr);
% goal. Find stim file with closest timing.
list=dir(stimDir);

stimFile='';
minDelta=100000;

% find closest stimFile in time to uMan start time (vStr)
for i=1:size(list,1)
    testStr=list(i).name;
    testVar=strsplit(testStr,'_');
    if str2num(testVar{1});
        numVar=str2num([testVar{1},testVar{2}]);
        testDelta=vStr-numVar;
        if abs(testDelta)<abs(minDelta)
            minDelta=testDelta;
            stimFile=testStr;
        end
    end
end

% Extract time & temp from stimFile

[ sTimeSec, sTemp ] = getStimInfo(stimFile,stimDir);

% get results from ROIs, using mean gray value ATM with an initial
% 1st measurement is background reading for subtraction...
% save without column headers
fVals= csvread(fullfile(imDir,'Results.csv'));
fVals=fVals(:,2:end);
BGnormVal=fVals-repmat(fVals(:,1),[1,size(fVals,2)]);
F0=repmat(min(fVals),[size(fVals,1),1]);
dFtoF=(fVals-F0)./F0;


% create time series sampled at 10hz to ease alignment
tsF=timeseries(dFtoF,vTimeSec);
tsT=timeseries(sTemp,sTimeSec);
[tsF,tsT] = synchronize(tsF,tsT,'Uniform','Interval',.1)



% Display options
diffdF=diff(dFtoF);
%set threshold... How?
figure(); histogram(diffdF);
threshDiff= 10*mean(median(abs(diffdF)));


for i=1:size(dFtoF,2)
    smDiff=diffdF(:,i); %smooth(diffdF(:,i),3);
    figure();hold on; plot(dFtoF(:,i)); yyaxis right; plot(smDiff);
    line(1:size(dFtoF,1),repmat(threshDiff,[size(dFtoF,1),1]));
end



% dFtoF x dT/dt
for i=1:size(dFtoF,2)
    smDiff=diff(tsT.data);
    figure();hold on; plot(tsF.data(:,i)); yyaxis right; plot(smDiff);
    line(1:size(dFtoF,1),repmat(threshDiff,[size(dFtoF,1),1]));
    text(.5,.9,strcat('worm #',num2str(i)),'units', 'normalized');
end

figure(); hold on;
for i=1:size(tsF.data,2)
    figure();plot(tsT.data,tsF.data(:,i),'linestyle', 'none','marker','*','color',rand(1,3))
    text(.5,.9,strcat('worm #',num2str(i)),'units', 'normalized');
end

figure(); hold on;
for i=1:size(tsF.data,2)
    figure(); plot(diff(tsT.data),diff(tsF.data(:,i)),'linestyle', 'none','marker','*','color',rand(1,3));
    
    text(.5,.9,strcat('worm #',num2str(i)),'units', 'normalized');
end