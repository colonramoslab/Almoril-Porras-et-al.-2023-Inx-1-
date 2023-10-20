% Tools & functions

% Pulling data into matlab
% need to write function, basically use CaImageAnalysis on all image files
% within a selected input parent directory, folderName
% or uiget parent directory if empty.
[outStatus] = caImageSubfolders(folderName)


% Alternatively, apply to single image stack
% Extracts time from video files, uiget driven.
% Aligns to fluorescence values in Results.csv file
% Creates & saves two synchronized time series, tsF & tsT, 
%       with values for Fluorescence & Temperature, respectively
% CaImageAnalysis_v0.m % replaced by function on 181120 to allow
% specification with uiget...

[tsF, tsT] = CaImageAnalysis('imageFile');

% Group combining tool
% pulls data from time series specified by dirs
% into single set of time series. 
% Fluorescence values and temperature values.
dirs={};
dirs{1}='C:\Users\jshha\Dropbox\ColonRamosLab\CaImaging_uMan\181009_DCR7450_1';
dirs{2}='C:\Users\jshha\Dropbox\ColonRamosLab\CaImaging_uMan\181009_DCR7450_2';
[temp_Wt, fluor_Wt] = loadTS(dirs);

% First pass figure, can include all frames or subset
incFrames=[710:820,1130:1240];
[ fM ] = makeHeatMapFigure( temp_Wt, fluor_Wt, 'incFrames', incFrames);

% Extract calcium transients, peaks
% optional screen to only analyze certain time periods.
% sampWin, amount of data to pull around identified peaks
% output is:
%       respMat: sample of fluorescence traces +/-sampWin around peaks
%       tempMat: sample of temperature traces +/-sampWin around peaks
%       respCnt: logical index into responses from fluor_Wt.data
[ screen] = makeScreen( fluor_Wt.data, [475:715, 2120:2360, 3745:3985] );
[ respMat, tempMat, respCnt] = respFind( fluor_Wt.data, temp_Wt.data,...
    'screen',screen,'sampWin',50);

% Response rate can be extracted by counting responses in respCnt
% divided by sampling interval time.
responseCount=sum(respCnt);
% here W1 is first bin of incFrames, ...
W1=[475:715]; W2=[2120:2360]; W3=[3745:3985]
sampDuration= fluor_Wt.Time(W1(end))- fluor_Wt.Time(W1(1))+...
    fluor_Wt.Time(W2(end))- fluor_Wt.Time(W2(1))+...
    fluor_Wt.Time(W3(end))- fluor_Wt.Time(W3(1));
RR=responseCount./sampDuration;

% This function may be better than piecemeal solutions:
incFrames={W1, W2, W3};
dirs=uigetdir(pwd,'directory with experiment for analysis:');
[ respMat, tempMat, respCnt, RR] = analyze_FreqAmp(dirs, incFrames)

% bring separate arrays together as NaN buffered matrix
% Makes it easier for figures.
[ dMat] = combineGroups({RR,RR});

% Figure type
% Makes a plot with each point from data set. Also circle at mean & error.
[fh] = plotEachPoint(dMat);

% Figure type
% Makes a single bar per group fragmented like pie chart with proportion in
% each group
 [fh] = barFractionPerGroup(dMat);

% Figure type
% Makes a plot with each point from data set. Also bar & error.
[fh] = plotEachPointBar(dMat)


% Saves figure as pdf in vector format for subsequent editing in
% illustrator, or similar vector graphics program
saveName='sampleFigure'; % can specify dir as well
 [ fStat, cnt ] = vectorSave( fh, saveName);

% Create figures comparing three groups within specified windows
winName='22to20C';
winCell={310:450,1520:1640,1940:2060,3150:3270,3560:3680}; 
sampNames={'Wt','e4-che2KO','AFDresc'};
[ fhs ] = compareGroupsWindow( fluor_Wt,fluor_che2KO,fluor_AFDresc,...
    temp_Wt,temp_che2KO,temp_AFDresc, winCell,...
    'sampNames',sampNames,'winName',winName);
 
 [ fhs ] = quantResponseConditions( G1F, G1T, C1,C2,C3)
 
% Finds responses, sPeaks, and extracts fluorescence & temperature information
% Useful for response-triggered averaging based on temperature.
% Could also look at AIY-response relationship to AFD flurescence or
% vice-versa
[ fluorMat, tempMat, sPeaks ] = respFind( fluor, temp, varargin )

% 

% %  % Example preliminary analysis
% %  Analysis_181022.m