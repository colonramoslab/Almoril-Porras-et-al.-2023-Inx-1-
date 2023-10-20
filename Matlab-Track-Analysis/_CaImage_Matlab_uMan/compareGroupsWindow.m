function [ fhs ] = compareGroupsWindow( G1F,G2F,G3F,G1T,G2T,G3T, winCell,varargin )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

winName='';
sampNames='';
winList=[]; winCnt=0;
G1time=0; G2time=0; G3time=0;

varargin=assignApplicable(varargin);

sampTimes=length(winCell);
for i=1:sampTimes
    tempWin=winCell{i};
    winList=[winList,tempWin];
    winCnt=winCnt+1;
    G1time=G1time+ G1T.time(tempWin(end))-G1T.time(tempWin(1));
    G2time=G2time+ G2T.time(tempWin(end))-G2T.time(tempWin(1));
    G3time=G3time+ G3T.time(tempWin(end))-G3T.time(tempWin(1));
end

[ screen] = makeScreen( G1F.data, winList);
[ respMat_G1, tempMat_G1, respCnt_G1] = respFind( G1F.data, G1T.data,'screen',screen);

[ screen] = makeScreen( G2F.data, winList);
[ respMat_G2, tempMat_G3, respCnt_G2] = respFind( G2F.data, G2T.data,'screen',screen);

[ screen] = makeScreen( G3F.data, winList);
[ respMat, tempMat, respCnt_G3] = respFind( G3F.data, G3T.data,'screen',screen);

% convert to Hz.
FR_G1=sum(respCnt_G1)./G1time;
FR_G2=sum(respCnt_G2)./G2time;
FR_G3=sum(respCnt_G3)./G3time;
[ dMat ] = combineGroups( {FR_G1,FR_G2,FR_G3} );

fhs{1} = barFractionPerGroup(dMat,'groupLabels',sampNames);
vectorSave(fhs{1},strcat('respBar_',winName))

fhs{2}= plotEachPoint(dMat,'groupLabels',sampNames);
vectorSave(fhs{2},strcat('FreqResp_',winName));


end

