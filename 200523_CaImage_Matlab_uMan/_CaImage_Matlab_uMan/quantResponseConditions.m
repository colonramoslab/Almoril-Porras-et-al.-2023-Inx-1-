function [ fhs ] = quantResponseConditions( G1F, G1T, C1,C2,C3,varargin)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

winName='';
sampNames='';
winList=[]; winCnt=0;
C1time=0; C2time=0; C3time=0;
C1list=[]; C2list=[]; C3list=[];

varargin=assignApplicable(varargin);

% Time in each phase
sampTimes=length(C1);
for i=1:sampTimes
    tempWin1=C1{i};
    C1list=[C1list,tempWin1];
    C1time=C1time+ G1T.time(tempWin1(end))-G1T.time(tempWin1(1));
end

sampTimes=length(C2);
for i=1:sampTimes
    tempWin2=C2{i};
    C2list=[C2list,tempWin2];
    C2time=C2time+ G1T.time(tempWin2(end))-G1T.time(tempWin2(1));
end

sampTimes=length(C3);
for i=1:sampTimes
    tempWin3=C3{i};
    C3list=[C3list,tempWin3];
    C3time=C3time+ G1T.time(tempWin3(end))-G1T.time(tempWin3(1));
end




% Count responses in each phase
[ screen] = makeScreen( G1F.data, C1list);
[ ~, ~, respCnt_C1] = respFind( G1F.data, G1T.data,'screen',screen);
[ screen] = makeScreen( G1F.data, C2list);
[ ~, ~, respCnt_C2] = respFind( G1F.data, G1T.data,'screen',screen);
[ screen] = makeScreen( G1F.data, C3list);
[ ~, ~, respCnt_C3] = respFind( G1F.data, G1T.data,'screen',screen);

% convert to Hz.
FR_G1=sum(respCnt_C1)./C1time;
FR_G2=sum(respCnt_C2)./C2time;
FR_G3=sum(respCnt_C3)./C3time;
[ dMat ] = combineGroups( {FR_G1,FR_G2,FR_G3} );


fhs{1} = barFractionPerGroup(dMat,'groupLabels',sampNames);
vectorSave(fhs{1},strcat('respBar_',winName))

fhs{2}= plotEachPoint(dMat,'groupLabels',sampNames);
vectorSave(fhs{2},strcat('FreqResp_',winName));


% Number of runs per condition type. Bar graph

% Fraction of runs. Pie chart

end

