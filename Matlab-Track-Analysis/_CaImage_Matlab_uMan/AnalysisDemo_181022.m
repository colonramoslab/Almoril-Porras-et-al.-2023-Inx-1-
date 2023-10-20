% Analysis_181022.m
% Analysis of caImaging data for eat-4 che-2p KO with and without
% AFD-specific rescue

% Make sure that the working directory is the parent of folder "demoData"

% % Bring relevant time series into workspace
PullSet_eat4KO % PullSet_eat4KO.m

saveDir=fullfile(pwd,'outputDir');
if ~exist(saveDir)
    mkdir(saveDir);
end

% %% visualize entire protocol
% incFrames=[1:1650];
% [ fM ] = makeHeatMapFigure( temp_Wt, fluor_Wt, 'incFrames', incFrames);
% vectorSave( fM, fullfile(saveDir,'Wt_first.pdf'));
% [ fM ] = makeHeatMapFigure( temp_che2KO, fluor_che2KO, 'incFrames', incFrames);
% vectorSave( fM, fullfile(saveDir,'che2p-e4KO_first.pdf'));
% [ fM ] = makeHeatMapFigure( temp_AFDresc, fluor_AFDresc, 'incFrames', incFrames);
% vectorSave( fM, fullfile(saveDir,'che2p-e4KO-AFDr_first.pdf'));


% %% first trial whole protocol
% 
% incFrames=[710:820,1130:1240];
% [ fM ] = makeHeatMapFigure( temp_Wt, fluor_Wt, 'incFrames', incFrames);
% vectorSave( fM, fullfile(saveDir,'Wt_rising1.pdf'));
% [ fM ] = makeHeatMapFigure( temp_che2KO, fluor_che2KO, 'incFrames', incFrames);
% vectorSave( fM, fullfile(saveDir,'che2p-e4KO_rising1.pdf'));
% [ fM ] = makeHeatMapFigure( temp_AFDresc, fluor_AFDresc, 'incFrames', incFrames);
% vectorSave( fM, fullfile(saveDir,'che2p-e4KO-AFDr_rising1.pdf'));
% 
% %% Narrow window on rising ramps, first trial only
% 
incFrames=[1520:1640,320:450];
[ fM ] = makeHeatMapFigure( temp_Wt, fluor_Wt, 'incFrames', incFrames);
vectorSave( fM, fullfile(saveDir,'Wt_cooling1.pdf'));
[ fM ] = makeHeatMapFigure( temp_che2KO, fluor_che2KO, 'incFrames', incFrames);
vectorSave( fM, fullfile(saveDir,'che2p-e4KO_cooling1.pdf'));
[ fM ] = makeHeatMapFigure( temp_AFDresc, fluor_AFDresc, 'incFrames', incFrames);
vectorSave( fM, fullfile(saveDir,'che2p-e4KO-AFDr_cooling1.pdf'));
 
%% To what are neurons responding?
% Response-locked temperature profile
% Find all responses. Align temperature before and after. 

% indicate responses, logical array 
% based on && of two features amplitude & derivative
% [ respMat, tempMat] = respFind( fluor_Wt.data, temp_Wt.data);
% % figure(); plot(mean(respMat,2))
% % figure(); plot(mean(tempMat,2))
% % figure(); plot(mean(diff(tempMat),2))

% [ respMat, tempMat ] = respFind( fluor_che2KO.data, temp_che2KO.data);
% % figure(); plot(mean(respMat,2))
% % figure(); plot(mean(tempMat,2))

% [ respMat, tempMat ] = respFind( fluor_AFDresc.data, temp_AFDresc.data);
% figure(); plot(mean(respMat,2))
% figure(); plot(mean(tempMat,2))

%% VERY clear that stimulus selection changes with che-2p eat-4KO & rescued
% in AFD

% % Is this true within plateaus? or just a consequence of ramp activity?
% mini-temperature change locked response profile?

% extract respMat & tempMat within certain intervals...
% make a matrix then use && overlay with respMat or tempMat




%% WT plateau response properties
%Plateau responses don't show a clear pattern. Undersampling or just ain't no pattern.

% wt during 18C holds
% [ screen] = makeScreen( fluor_Wt.data, [475:715, 2120:2360, 3745:3985] );
% [ respMat, tempMat] = respFind( fluor_Wt.data, temp_Wt.data,'screen',screen);
% figure(); plot(mean(respMat,2))
% figure(); plot(mean(tempMat,2))
% figure(); plot(mean(diff(tempMat),2))

% wt during 20C holds
%     scWin_20C_Initial=[50:290];
%     scWin_20C_postLow=[880:1120,2500:2740, 3320:3560];
%     scWin_20C_postHigh=[1700:1940, 3320:3560, 4140:4380];
% [ screen] = makeScreen( fluor_Wt.data,[scWin_20C_Initial,scWin_20C_postLow,scWin_20C_postHigh] );
% [ respMat, tempMat] = respFind( fluor_Wt.data, temp_Wt.data,'screen',screen);
% figure(); plot(mean(respMat,2))
% figure(); plot(mean(tempMat,2))
% figure(); plot(mean(diff(tempMat),2))

% wt during 22C hold
% [ screen] = makeScreen( fluor_Wt.data,[1270:1510, 2910:3150, 4560:4800]);
% [ respMat, tempMat] = respFind( fluor_Wt.data, temp_Wt.data,'screen',screen);
% figure(); plot(mean(respMat,2))
% figure(); plot(mean(tempMat,2))
% figure(); plot(mean(diff(tempMat),2))


%% eat-4 plateau response properties
% che-2 eat-4 KO during 18C holds
% [ screen] = makeScreen( fluor_che2KO.data, scWin_18C );
% [ respMat, tempMat] = respFind( fluor_che2KO.data, temp_che2KO.data,'screen',screen);
% figure(); plot(mean(respMat,2))
% figure(); plot(mean(tempMat,2))
% figure(); plot(mean(diff(tempMat),2))

% che-2 during 20C holds
% [ screen] = makeScreen( fluor_che2KO.data,[50:290,880:1120,2500:2740, 3300:3540,1680:1920, 3305:3545, 4115:4345] );
% [ respMat, tempMat] = respFind( fluor_che2KO.data, temp_che2KO.data,'screen',screen);
% figure(); plot(mean(respMat,2))
% figure(); plot(mean(tempMat,2))
% figure(); plot(tempMat)
% figure(); plot(mean(diff(tempMat),2))

% che-2 eat-4 KO during 22C hold
% [ screen] = makeScreen( fluor_che2KO.data,[1270:1510, 2910:3150, 4560:4800]);
% [ respMat, tempMat] = respFind( fluor_che2KO.data, temp_che2KO.data,'screen',screen);
% figure(); plot(mean(respMat,2))
% figure(); plot(mean(tempMat,2))
% figure(); plot(tempMat)
% figure(); plot(mean(diff(tempMat),2))

%% eat-4 che-2 KO + AFD::eat-4 plateau response properties
% che-2 eat-4 KO with AFD rescue during 18C holds
% [ screen] = makeScreen( fluor_AFDresc.data, [475:715, 2120:2360,3715:3955]); 
% [ respMat, tempMat] = respFind( fluor_AFDresc.data, temp_AFDresc.data,'screen',screen);
% figure(); plot(tempMat)
% figure(); plot(mean(respMat,2))
% figure(); plot(mean(tempMat,2))
% figure(); plot(mean(diff(tempMat),2))

% che-2 with AFD rescue during 20C holds
% [ screen] = makeScreen( fluor_AFDresc.data,[70:290, 910:1150,2560:2800,3300:3540, 1680:1920,3305:3545, 4200:4420] );
% [ respMat, tempMat] = respFind( fluor_AFDresc.data, temp_AFDresc.data,'screen',screen);
% figure(); plot(tempMat)
% figure(); plot(mean(respMat,2))
% figure(); plot(mean(tempMat,2))
% figure(); plot(mean(diff(tempMat),2))

% che-2 eat-4 KO with AFD rescue during 22C hold
% [ screen] = makeScreen( fluor_AFDresc.data,[1270:1510, 2960:3140, 4560:4770]);
% [ respMat, tempMat] = respFind( fluor_AFDresc.data, temp_AFDresc.data,'screen',screen);
% figure(); plot(tempMat)
% figure(); plot(mean(respMat,2))
% figure(); plot(mean(tempMat,2))
% figure(); plot(mean(diff(tempMat),2))

%% What about warming ramps?

% Do the same for ramps. At least this will give a response count for
% ramping intervals.

W1=[715:835];
W2=[2350:2470];
W3=[3970:4090];
% wt 18C->20C. 3 responses from 
[ screen] = makeScreen( fluor_Wt.data, [W1,W2,W3]);
[ respMat, tempMat, respCnt_wt] = respFind( fluor_Wt.data, temp_Wt.data,'screen',screen);
% figure(); plot(tempMat)
% figure(); plot(mean(respMat,2))
% figure(); plot(mean(tempMat,2))
% figure(); plot(mean(diff(tempMat),2))

% eat-4 che-2pKO 18C->20C. 14 responses from 16worm over 3 trials
[ screen] = makeScreen( fluor_che2KO.data,  [W1,W2,W3]);
[ respMat, tempMat, respCnt_che2KO] = respFind( fluor_che2KO.data, temp_che2KO.data,'screen',screen);
% figure(); plot(tempMat)
% figure(); plot(mean(respMat,2))
% figure(); plot(mean(tempMat,2))
% figure(); plot(mean(diff(tempMat),2))

% eat-4 che-2pKO + AFDeat4 18C->20C. 6 responses from 
[ screen] = makeScreen( fluor_AFDresc.data,  [W1,W2,W3]);
[ respMat, tempMat, respCnt_AFDresc] = respFind( fluor_AFDresc.data, temp_AFDresc.data,'screen',screen);
% figure(); plot(tempMat)
% figure(); plot(respMat)
% figure(); plot(mean(respMat,2))
% figure(); plot(mean(tempMat,2))
% figure(); plot(mean(diff(tempMat),2))

RR_wt=sum(respCnt_wt)./3;
RR_che2KO=sum(respCnt_che2KO)./3;
RR_AFDresc=sum(respCnt_AFDresc)./3;

[ dMat ] = combineGroups( {RR_wt,RR_che2KO,RR_AFDresc} );
fh= plotEachPoint(dMat);
[ fh ] = barFractionPerGroup(dMat);

% Should all be the same sampling time, but confirm.
sampInt_wt= fluor_Wt.Time(W1(end))- fluor_Wt.Time(W1(1))+...
    fluor_Wt.Time(W2(end))- fluor_Wt.Time(W2(1))+...
    fluor_Wt.Time(W3(end))- fluor_Wt.Time(W3(1));

sampInt_che2KO= fluor_che2KO.Time(W1(end))- fluor_che2KO.Time(W1(1))+...
    fluor_che2KO.Time(W2(end))- fluor_che2KO.Time(W2(1))+...
    fluor_che2KO.Time(W3(end))- fluor_che2KO.Time(W3(1));

sampInt_AFDresc= fluor_AFDresc.Time(W1(end))- fluor_AFDresc.Time(W1(1))+...
    fluor_AFDresc.Time(W2(end))- fluor_AFDresc.Time(W2(1))+...
    fluor_AFDresc.Time(W3(end))- fluor_AFDresc.Time(W3(1));

% convert to Hz.
FR_wt=sum(respCnt_wt)./sampInt_wt;
FR_che2KO=sum(respCnt_che2KO)./sampInt_che2KO;
FR_AFDresc=sum(respCnt_AFDresc)./sampInt_AFDresc;
[ dMat ] = combineGroups( {FR_wt,FR_che2KO,FR_AFDresc} );
fh= plotEachPoint(dMat,'FreqResp_18to20C',{'Wt','eat4-Che2KO','AFDresc'});

%% Second warming ramp
% % To define window, plot the temperature profile for all.
% % figure(); hold on; plot(temp_Wt.data); plot(temp_che2KO.data); plot(temp_AFDresc.data);
winName='20to22C';
winCell={1110:1240, 2740:2870, 4390:4520};
sampNames={'Wt','e4-che2KO','AFDresc'};
[ fhs ] = compareGroupsWindow( fluor_Wt,fluor_che2KO,fluor_AFDresc,...
    temp_Wt,temp_che2KO,temp_AFDresc, winCell,...
    'sampNames',sampNames,'winName',winName);

figure(); hold on; plot(temp_Wt.data); plot(temp_che2KO.data); plot(temp_AFDresc.data);

%% BOTH warming ramps
% % To define window, plot the temperature profile for all.
% % figure(); hold on; plot(temp_Wt.data); plot(temp_che2KO.data); plot(temp_AFDresc.data);
winName='BothWarming';
winCell={715:835,1110:1240,2350:2470, 2740:2870, 3970:4090,4390:4520};
sampNames={'Wt','e4-che2KO','AFDresc'};
[ fhs ] = compareGroupsWindow( fluor_Wt,fluor_che2KO,fluor_AFDresc,...
    temp_Wt,temp_che2KO,temp_AFDresc, winCell,...
    'sampNames',sampNames,'winName',winName);

figure(); hold on; plot(temp_Wt.data); plot(temp_che2KO.data); plot(temp_AFDresc.data);


%% 18C plateaus
% % To define window, plot the temperature profile for all.
% % figure(); hold on; plot(temp_Wt.data); plot(temp_che2KO.data); plot(temp_AFDresc.data);
winName='18Cplat';
winCell={475:715, 2120:2360, 3745:3985} ;
sampNames={'Wt','e4-che2KO','AFDresc'};
[ fhs ] = compareGroupsWindow( fluor_Wt,fluor_che2KO,fluor_AFDresc,...
    temp_Wt,temp_che2KO,temp_AFDresc, winCell,...
    'sampNames',sampNames,'winName',winName);

%% 20C plateaus
% % To define window, plot the temperature profile for all.
% % figure(); hold on; plot(temp_Wt.data); plot(temp_che2KO.data); plot(temp_AFDresc.data);
winName='20Cplat';
winCell={50:290,880:1120,2500:2740, 3320:3560,1700:1940, 3320:3560, 4140:4380};
sampNames={'Wt','e4-che2KO','AFDresc'};
[ fhs ] = compareGroupsWindow( fluor_Wt,fluor_che2KO,fluor_AFDresc,...
    temp_Wt,temp_che2KO,temp_AFDresc, winCell,...
    'sampNames',sampNames,'winName',winName);

%% 22C plateaus
% % To define window, plot the temperature profile for all.
% % figure(); hold on; plot(temp_Wt.data); plot(temp_che2KO.data); plot(temp_AFDresc.data);
winName='22Cplat';
winCell={1270:1510, 2960:3140, 4560:4770};
sampNames={'Wt','e4-che2KO','AFDresc'};
[ fhs ] = compareGroupsWindow( fluor_Wt,fluor_che2KO,fluor_AFDresc,...
    temp_Wt,temp_che2KO,temp_AFDresc, winCell,...
    'sampNames',sampNames,'winName',winName);

%% Cooling ramps, 20to18C
% % To define window, plot the temperature profile for all.
% % figure(); hold on; plot(temp_Wt.data); plot(temp_che2KO.data); plot(temp_AFDresc.data);
winName='20to18C';
winCell={310:450, 1920:2060, 3540:3680};
sampNames={'Wt','e4-che2KO','AFDresc'};
[ fhs ] = compareGroupsWindow( fluor_Wt,fluor_che2KO,fluor_AFDresc,...
    temp_Wt,temp_che2KO,temp_AFDresc, winCell,...
    'sampNames',sampNames,'winName',winName);


%% Cooling ramps, 22to20C
% % To define window, plot the temperature profile for all.
% % figure(); hold on; plot(temp_Wt.data); plot(temp_che2KO.data); plot(temp_AFDresc.data);
winName='22to20C';
winCell={1520:1640, 3150:3270};%, 4760:4850};
sampNames={'Wt','e4-che2KO','AFDresc'};
[ fhs ] = compareGroupsWindow( fluor_Wt,fluor_che2KO,fluor_AFDresc,...
    temp_Wt,temp_che2KO,temp_AFDresc, winCell,...
    'sampNames',sampNames,'winName',winName);

%% Both Cooling ramps
% % To define window, plot the temperature profile for all.
% % figure(); hold on; plot(temp_Wt.data); plot(temp_che2KO.data); plot(temp_AFDresc.data);
winName='22to20C';
winCell={310:450,1520:1640,1940:2060,3150:3270,3560:3680}; 
sampNames={'Wt','e4-che2KO','AFDresc'};
[ fhs ] = compareGroupsWindow( fluor_Wt,fluor_che2KO,fluor_AFDresc,...
    temp_Wt,temp_che2KO,temp_AFDresc, winCell,...
    'sampNames',sampNames,'winName',winName);


%% Compare within condition between sample periods, wild-type
warmC={715:835,1110:1240,2350:2470, 2740:2870, 3970:4090,4390:4520};
platC={50:290,880:1120,2500:2740, 3320:3560,1700:1940, 3320:3560, 4140:4380};
coolC={310:450,1520:1640,1940:2060,3150:3270,3560:3680};
winName='Wt';
sampNames={'Warm','Plat','Cool'};
[ fhs ] = quantResponseConditions( fluor_Wt, temp_Wt, warmC,platC,coolC,...
    'sampNames',sampNames,'winName',winName);

%% Compare within condition between sample periods, eat-4 che-2KO
warmC={715:835,1110:1240,2350:2470, 2740:2870, 3970:4090,4390:4520};
platC={50:290,880:1120,2500:2740, 3320:3560,1700:1940, 3320:3560, 4140:4380};
coolC={310:450,1520:1640,1940:2060,3150:3270,3560:3680};
winName='che2pKO';
sampNames={'Warm','Plat','Cool'};
[ fhs ] = quantResponseConditions( fluor_che2KO, temp_che2KO, warmC,platC,coolC,...
    'sampNames',sampNames,'winName',winName);


%% Compare within condition between sample periods, eat-4 che-2KO + AFDeat4
warmC={715:835,1110:1240,2350:2470, 2740:2870, 3970:4090,4390:4520};
platC={50:290,880:1120,2500:2740, 3320:3560,1700:1940, 3320:3560, 4140:4380};
coolC={310:450,1520:1640,1940:2060,3150:3270,3560:3680};
winName='AFDe4resc';
sampNames={'Warm','Plat','Cool'};
[ fhs ] = quantResponseConditions( fluor_AFDresc, temp_AFDresc, warmC,platC,coolC,...
    'sampNames',sampNames,'winName',winName);